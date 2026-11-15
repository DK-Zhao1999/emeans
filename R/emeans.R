#' Batch size and step size from the e-means scaling law
#'
#' The paper fixes the batch size and the initial step size as a function of
#' the problem shape, so that the number of updates stays roughly constant as
#' the data grow:
#'
#' \deqn{b = c_b \, k^{0.13} \, n^{0.88} \, E^{0.81}}
#' \deqn{\alpha_0 = 0.824 \, k^{-0.46} \, b^{0.90} \, n^{-0.45} \, E^{-0.11}}
#'
#' The batch size is rounded to a power of two and clamped to `[16, n]`; the
#' step size is capped at 1.
#'
#' @param n Number of samples.
#' @param k Number of clusters.
#' @param epochs Number of passes over the data, default 5.
#' @param quality Logical. Use the higher-quality constant
#'   (\eqn{c_b = 0.0044}, far more, smaller batches) instead of the default
#'   \eqn{c_b = 0.035}.
#'
#' @return A list with elements
#'   * `batch`: integer batch size, a power of two in `[16, n]`;
#'   * `alpha0`: initial step size, in `(0, 1]`;
#'   * `capped`: whether the law asked for a step size above 1 and it was
#'     clipped.
#'
#' @examples
#' emeans_law(n = 500, k = 3)
#' emeans_law(n = 5000, k = 2, quality = TRUE)$batch
#'
#' @export
emeans_law <- function(n, k, epochs = 5, quality = FALSE) {
  c_b <- if (quality) 0.0044 else 0.035
  b <- c_b * k^0.13 * n^0.88 * epochs^0.81
  if (n > 2^28) b <- b / 8
  b <- 2^round(log2(b))
  b <- min(max(16, b), n)                    # clamp to [16, n]: max first, then min

  a0 <- 0.824 * k^(-0.46) * b^0.90 * n^(-0.45) * epochs^(-0.11)
  capped <- a0 > 1
  a0 <- min(a0, 1)

  list(batch = as.integer(b), alpha0 = a0, capped = capped)
}


## -----------------------------------------------------------------------------
## Nearest centroid, in row chunks so the distance matrix stays small
## -----------------------------------------------------------------------------
.emeans_nearest <- function(xb, C, C2 = NULL, chunk = 8192L) {
  if (is.null(C2)) C2 <- rowSums(C * C)
  nb <- nrow(xb)
  labels <- integer(nb)
  dmin <- numeric(nb)
  for (s in seq.int(1L, nb, by = chunk)) {
    e <- min(nb, s + chunk - 1L)
    xs <- xb[s:e, , drop = FALSE]
    d2 <- -2 * tcrossprod(xs, C) + rowSums(xs * xs)
    d2 <- d2 + rep(C2, each = e - s + 1L)
    lab <- max.col(-d2, ties.method = "first")
    labels[s:e] <- lab
    dmin[s:e] <- d2[cbind(seq_len(e - s + 1L), lab)]
  }
  list(labels = labels, dmin = dmin)
}


#' e-means clustering engine
#'
#' Clusters rows of a numeric matrix with e-means. This is the low-level
#' engine: it expects the matrix to be already preprocessed (samples in rows,
#' features in columns) and receives `k` explicitly. For expression matrices
#' use [emeans_bio()], which does the preprocessing, the choice of `k` and the
#' mapping back to genes.
#'
#' @param x Numeric matrix, one sample per row and one feature per column.
#' @param k Number of clusters, a positive integer. Values above `nrow(x)` are
#'   reduced to `nrow(x)` with a warning.
#' @param epochs Number of passes over the data, default 5. Two or three are
#'   usually enough for the centroids to settle.
#' @param batch Batch size. `NULL` (default) takes it from [emeans_law()];
#'   a user value is clamped to `[1, nrow(x)]`.
#' @param alpha0 Initial step size in `(0, 1]`. `NULL` (default) takes it from
#'   [emeans_law()].
#' @param schedule `"cosine"` (default) anneals the step size from `alpha0`
#'   down to `rho * alpha0`; `"constant"` keeps it at `alpha0`.
#' @param rho Floor of the cosine schedule, as a fraction of `alpha0`.
#' @param init Optional `k` by `ncol(x)` matrix of starting centroids. `NULL`
#'   (default) samples `k` rows of `x`.
#' @param revival Logical; split under-used centroids. Strongly recommended.
#' @param revival_every Revival is attempted every this many updates.
#' @param revival_frac A centroid whose weight falls below this fraction of the
#'   mean weight is considered dead.
#' @param finalize How to produce the returned centroids: `"last-epoch-mean"`
#'   (default) recomputes them as the exact mean of the samples assigned in the
#'   last epoch, `"slow-ema"` uses a long-running average, `"none"` keeps the
#'   final EMA state.
#' @param slow_half_life Half-life, in updates, of the slow average used by
#'   `finalize = "slow-ema"`. Defaults to a quarter of the total updates.
#' @param seed Optional random seed.
#' @param verbose Print progress.
#'
#' @return An object of class `emeans`, a list with elements `centroids`
#'   (`k` by `ncol(x)`), `sizes`, `labels`, `mse`, `batch`, `alpha0`,
#'   `updates`, `epochs`, `k`, `n`, `revived`, `elapsed` and `history`
#'   (a data frame with one row per update).
#'
#' @examples
#' set.seed(1)
#' x <- rbind(matrix(rnorm(200, 0), 100, 2), matrix(rnorm(200, 6), 100, 2))
#' fit <- emeans(x, k = 2, seed = 1)
#' fit
#' table(fit$labels)
#'
#' @seealso [emeans_bio()] for the expression-matrix interface.
#'
#' @export
emeans <- function(x, k, epochs = 5L,
                   batch = NULL, alpha0 = NULL,
                   schedule = c("cosine", "constant"),
                   rho = 0.1,
                   init = NULL,
                   revival = TRUE, revival_every = 5L, revival_frac = 0.1,
                   finalize = c("last-epoch-mean", "slow-ema", "none"),
                   slow_half_life = NULL,
                   seed = NULL, verbose = FALSE) {

  schedule <- match.arg(schedule)
  finalize <- match.arg(finalize)
  if (!is.null(seed)) set.seed(seed)

  x <- as.matrix(x)
  if (!is.numeric(x)) stop("`x` must be a numeric matrix", call. = FALSE)
  storage.mode(x) <- "double"
  if (anyNA(x)) stop("`x` must not contain NA; use `emeans_bio()` to impute",
                     call. = FALSE)
  n <- nrow(x); d <- ncol(x)
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 1L) {
    stop("`k` must be a positive integer", call. = FALSE)
  }
  if (k > n) {
    warning("k > n: k is set to n", call. = FALSE)
    k <- n
  }
  epochs <- as.integer(epochs)
  if (length(epochs) != 1L || is.na(epochs) || epochs < 1L) {
    stop("`epochs` must be >= 1", call. = FALSE)
  }

  law <- emeans_law(n, k, epochs)
  law_alpha <- is.null(alpha0)
  law_batch <- is.null(batch)
  if (law_batch) batch <- law$batch
  ## a user-supplied batch only has to be in [1, n]; the law's floor of 16
  ## applies when the value comes from the law itself
  batch <- as.integer(if (law_batch) min(max(16L, batch), n) else min(max(1L, batch), n))
  if (law_alpha) alpha0 <- law$alpha0
  if (!is.finite(alpha0) || alpha0 <= 0 || alpha0 > 1) {
    stop("`alpha0` must be in (0, 1]", call. = FALSE)
  }

  n_batch <- ceiling(n / batch)
  total_updates <- epochs * n_batch
  if (is.null(slow_half_life)) slow_half_life <- max(2, total_updates / 4)
  alpha_slow <- 1 - exp(-log(2) / slow_half_life)

  alpha_at <- function(t) {
    if (schedule == "constant") return(alpha0)
    if (total_updates <= 1L) return(alpha0)
    alpha0 * (rho + (1 - rho) / 2 * (1 + cos(pi * t / (total_updates - 1))))
  }

  if (is.null(init)) {
    C <- x[sample.int(n, k, replace = FALSE), , drop = FALSE]
  } else {
    C <- as.matrix(init)
    if (!identical(dim(C), c(k, d))) {
      stop("`init` must be a k x d matrix", call. = FALSE)
    }
    storage.mode(C) <- "double"
  }
  ## start with N = 1e-10 so that S / N reproduces the given centroids exactly
  N <- rep(1e-10, k)
  S <- C * N
  C2 <- rowSums(C * C)
  Ns <- N; Ss <- S

  if (verbose) {
    cat(sprintf("e-means: n=%d d=%d k=%d epochs=%d batch=%d updates=%d alpha0=%.4g%s%s\n",
                n, d, k, epochs, batch, total_updates, alpha0,
                if (isTRUE(law$capped) && law_alpha) " (law capped at 1)" else "",
                if (!law_batch) " (user batch)" else ""))
  }

  hist_upd <- integer(total_updates); hist_ep <- integer(total_updates)
  hist_a <- numeric(total_updates); hist_mse <- numeric(total_updates)
  h <- 0L
  last_lab <- integer(n)
  n_revived <- 0L
  t0 <- proc.time()[["elapsed"]]

  t_idx <- 0L
  for (ep in seq_len(epochs)) {
    perm <- sample.int(n)
    for (b in seq_len(n_batch)) {
      sel <- perm[((b - 1L) * batch + 1L):min(n, b * batch)]
      xb <- x[sel, , drop = FALSE]

      nn <- .emeans_nearest(xb, C, C2)

      n_j <- tabulate(nn$labels, k)
      s_j <- matrix(0, k, d)
      rs <- rowsum(xb, group = nn$labels, reorder = TRUE)
      g <- as.integer(rownames(rs))
      s_j[g, ] <- rs

      a <- alpha_at(t_idx)
      t_idx <- t_idx + 1L
      N <- (1 - a) * N + a * n_j
      S <- (1 - a) * S + a * s_j
      ## when alpha0 is 1 and a cluster gets no sample, N drops to exactly 0
      ## and the centroid would silently collapse to the origin
      C <- S / pmax(N, 1e-30)
      bad <- which(!is.finite(rowSums(C)) | N <= 0)
      if (length(bad)) {
        C[bad, ] <- xb[sample.int(nrow(xb), length(bad), replace = TRUE), , drop = FALSE]
        N[bad] <- 1e-10
        S[bad, ] <- C[bad, ] * 1e-10
      }
      C2 <- rowSums(C * C)

      Ns <- (1 - alpha_slow) * Ns + alpha_slow * n_j
      Ss <- (1 - alpha_slow) * Ss + alpha_slow * s_j

      ## revival: every `revival_every` updates, split centroids whose weight
      ## fell below `revival_frac` of the mean
      if (isTRUE(revival) && (t_idx %% revival_every == 0L)) {
        thr <- revival_frac * mean(N)
        dead <- which(N < thr)
        live <- setdiff(seq_len(k), dead)
        if (length(dead) && length(live)) {
          ## sample.int() rather than sample(live, ...): when `live` has length
          ## 1, sample() would read it as a size and mis-match `prob`
          pick <- sample.int(length(live), size = length(dead), replace = TRUE,
                             prob = N[live])
          donor <- live[pick]
          n_revived <- n_revived + length(dead)
          for (ii in seq_along(dead)) {
            j <- dead[ii]; dd <- donor[ii]
            p <- xb[sample.int(nrow(xb), 1L), ]
            delta <- 0.1 * (p - C[dd, ])
            if (all(delta == 0)) delta <- rnorm(d, sd = 1e-3)
            c_don <- C[dd, ]
            half <- N[dd] / 2
            C[dd, ] <- c_don - delta; N[dd] <- half; S[dd, ] <- (c_don - delta) * half
            C[j,  ] <- c_don + delta; N[j]  <- half; S[j,  ] <- (c_don + delta) * half
          }
          C2 <- rowSums(C * C)
        }
      }

      h <- h + 1L
      hist_upd[h] <- t_idx
      hist_ep[h] <- ep
      hist_a[h] <- a
      hist_mse[h] <- mean(nn$dmin)

      if (ep == epochs) last_lab[sel] <- nn$labels
    }
    if (verbose && epochs > 1) {
      cat(sprintf("  epoch %d/%d  elapsed %.1fs  batch-mse %.6g\n",
                  ep, epochs, proc.time()[["elapsed"]] - t0,
                  mean(hist_mse[hist_ep == ep])))
    }
  }

  ## ---- finalization ----
  sizes <- tabulate(last_lab, k)
  if (finalize == "last-epoch-mean") {
    C_out <- C
    if (sum(sizes > 0)) {
      rs <- rowsum(x, group = last_lab, reorder = TRUE)
      g <- as.integer(rownames(rs))
      C_out[g, ] <- rs / sizes[g]
    }
  } else if (finalize == "slow-ema") {
    C_out <- Ss / pmax(Ns, 1e-30)
    bad <- !is.finite(rowSums(C_out)) | sizes == 0
    if (any(bad)) C_out[bad, ] <- C[bad, ]
  } else {
    C_out <- C
  }

  fit <- .emeans_nearest(x, C_out)

  structure(list(
    centroids = C_out,
    sizes     = tabulate(fit$labels, k),
    labels    = fit$labels,
    mse       = mean(fit$dmin),
    batch     = batch,
    alpha0    = alpha0,
    updates   = total_updates,
    epochs    = epochs,
    k         = k,
    n         = n,
    revived   = n_revived,
    elapsed   = proc.time()[["elapsed"]] - t0,
    history   = data.frame(update = hist_upd, epoch = hist_ep,
                           alpha = hist_a, batch_mse = hist_mse)
  ), class = "emeans")
}


#' @param x An `emeans` fit.
#' @param ... Ignored.
#' @export
#' @rdname emeans
print.emeans <- function(x, ...) {
  cat(sprintf("<e-means fit> n=%d k=%d epochs=%d batch=%d alpha0=%.3g updates=%d\n",
              x$n, x$k, x$epochs, x$batch, x$alpha0, x$updates))
  cat(sprintf("  train MSE = %.6g  elapsed = %.2fs  dead clusters = %d\n",
              x$mse, x$elapsed, sum(x$sizes == 0)))
  invisible(x)
}


#' Assign rows of a matrix to existing centroids
#'
#' @param x Numeric matrix, samples in rows.
#' @param centroids Numeric matrix, one centroid per row, same columns as `x`.
#'
#' @return A list with `labels`, the mean squared distance `mse` and the
#'   per-sample squared distance to its centroid (`dist`).
#'
#' @examples
#' x <- rbind(matrix(rnorm(40, 0), 20, 2), matrix(rnorm(40, 5), 20, 2))
#' cen <- rbind(c(0, 0), c(5, 5))
#' emeans_assign(x, cen)$labels
#'
#' @export
emeans_assign <- function(x, centroids) {
  x <- as.matrix(x); storage.mode(x) <- "double"
  C <- as.matrix(centroids); storage.mode(C) <- "double"
  if (ncol(x) != ncol(C)) {
    stop("`x` and `centroids` must have the same number of columns", call. = FALSE)
  }
  nn <- .emeans_nearest(x, C)
  list(labels = nn$labels, mse = mean(nn$dmin), dist = nn$dmin)
}
