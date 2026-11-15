#' Unsupervised subtype discovery for expression data
#'
#' Takes a genes-by-samples expression matrix, preprocesses it, clusters the
#' samples with [emeans()] and returns per-sample subtype labels together with
#' expression centroids, diagnostics and the pieces needed to classify new
#' samples.
#'
#' @details
#' The pipeline, in order:
#'
#' 1. samples are moved to rows and genes to columns;
#' 2. missing values are imputed with the gene median, genes that are missing
#'    in every sample are dropped;
#' 3. `Inf`/`-Inf` trigger an error rather than being removed quietly;
#' 4. zero-variance genes are dropped;
#' 5. each gene is centred and scaled (z-score);
#' 6. the first `pca` principal components are kept (noise reduction);
#' 7. [emeans()] is run on the scores, `n_init` random starts per `k`;
#' 8. `k` is chosen by mean silhouette width, unless it is given;
#' 9. centroids are recomputed in the gene space as exact subtype means.
#'
#' Clustering happens in the PC space, so the reported `silhouette` is measured
#' there as well and is not numerically comparable with a silhouette computed
#' on the raw expression distances.
#'
#' @param x Numeric matrix or data frame, **genes in rows and samples in
#'   columns** (the usual expression layout). Row and column names, if present,
#'   are carried through.
#' @param k Number of subtypes. `NULL` (default) selects it from `2:k_max` by
#'   mean silhouette width. A single integer of at least 2 fixes it.
#' @param k_max Upper limit of the automatic search, default 6.
#' @param epochs Passes over the data per run, default 5; 2-3 is usually enough.
#' @param n_init Number of random starts per candidate `k`; the run with the
#'   lowest training error is kept. Default 10.
#' @param pca Number of principal components to cluster on, default 50.
#'   `NULL` clusters on all genes.
#' @param center,scale Logical; centre and z-score each gene. Default both
#'   `TRUE`.
#' @param finalize How the final centroids are produced, see [emeans()].
#'   `"slow-ema"` tends to generalize better to new samples.
#' @param seed Random seed, default 1. `NULL` also means 1, so the result is
#'   reproducible by default.
#' @param verbose Print progress, default `TRUE`.
#'
#' @return An object of class `emeans_bio`, a list with elements
#'   * `cluster`: factor of subtype labels, named with the sample names;
#'   * `k`, `sizes`: number of subtypes and samples per subtype;
#'   * `mse`: training mean squared distance in the PC space;
#'   * `silhouette`: mean silhouette width, `silhouette_values` per (possibly
#'     subsampled) sample, `silhouette_idx` the samples used;
#'   * `k_scan`: data frame with `k`, `silhouette`, `mse`, `seconds`, one row
#'     per candidate `k`;
#'   * `centroids`: subtype by gene matrix of means, in the scaled space;
#'   * `centroids_expr`: the same centroids back on the input scale, which is
#'     what you want for biological interpretation;
#'   * `centroids_pc`: subtype by PC centroids, used by `predict()`;
#'   * `scores`: sample by PC coordinates;
#'   * `pve`: proportion of variance explained by each retained PC;
#'   * `fit`: the underlying `emeans` object;
#'   * `prep`: preprocessing parameters (gene names, centring, scaling, the
#'     rotation, the training medians);
#'   * `expr_scaled`: sample by gene scaled matrix;
#'   * `n_init`, `call`, `elapsed`.
#'
#' @examples
#' fit <- emeans_bio(emeans_example, seed = 1, verbose = FALSE)
#' fit
#' table(fit$cluster)
#' head(emeans_markers(fit, n = 5))
#'
#' @seealso [emeans_markers()], [predict.emeans_bio()], [plot.emeans_bio()].
#'
#' @export
emeans_bio <- function(x, k = NULL, k_max = 6L, epochs = 5L, n_init = 10L,
                       pca = 50, center = TRUE, scale = TRUE,
                       finalize = c("last-epoch-mean", "slow-ema"),
                       seed = 1L, verbose = TRUE) {
  finalize <- match.arg(finalize)
  if (is.null(seed)) seed <- 1L
  t_start <- proc.time()[["elapsed"]]

  prep <- .emeans_prep(x, center = center, scale = scale, pca = pca)
  Z <- prep$scores
  n <- nrow(Z)

  if (is.null(k)) {
    k_max <- max(2L, min(as.integer(k_max), n - 1L))
    ks <- 2L:k_max
  } else {
    ks <- as.integer(k)
    if (length(ks) != 1L || is.na(ks) || ks < 2L) {
      stop("`k` must be a single integer >= 2 (or NULL for automatic selection)",
           call. = FALSE)
    }
    if (ks > n - 1L) stop("`k` must be < number of samples", call. = FALSE)
  }
  if (verbose) {
    cat(sprintf("e-means subtyping: %d samples x %d genes -> %s PCs\n",
                n, length(prep$gene),
                if (is.na(prep$npc)) "no" else as.character(prep$npc)))
    if (is.null(k)) {
      cat(sprintf("  scanning k = %s (silhouette), %d inits each\n",
                  paste(ks, collapse = "-"), n_init))
    }
  }

  ## ---- run n_init starts for every k and keep the lowest training error ----
  best <- vector("list", length(ks))
  scan <- data.frame(k = integer(0), silhouette = numeric(0),
                     mse = numeric(0), seconds = numeric(0))
  for (i in seq_along(ks)) {
    kk <- ks[i]
    t0 <- proc.time()[["elapsed"]]
    cand <- NULL
    for (r in seq_len(max(1L, n_init))) {
      ## keep seed + r inside the integer range
      s_r <- as.integer((as.numeric(seed) + r - 1) %% .Machine$integer.max)
      f <- emeans(Z, k = kk, epochs = epochs, finalize = finalize,
                  seed = s_r, verbose = FALSE)
      if (is.null(cand) || f$mse < cand$mse) cand <- f
    }
    sil <- .emeans_silhouette(Z, cand$labels)
    best[[i]] <- cand
    scan <- rbind(scan, data.frame(k = kk, silhouette = sil$mean, mse = cand$mse,
                                   seconds = proc.time()[["elapsed"]] - t0))
    if (verbose && length(ks) > 1) {
      cat(sprintf("  k=%d  silhouette %.3f  train MSE %.4g  (%.1fs)\n",
                  kk, sil$mean, cand$mse, scan$seconds[i]))
    }
  }

  ## ---- choose k ----
  if (is.null(k)) {
    ok <- is.finite(scan$silhouette)
    if (!any(ok)) {
      sel <- 1L
      warning("silhouette undefined, falling back to smallest k", call. = FALSE)
    } else {
      sel <- which.max(ifelse(ok, scan$silhouette, -Inf))
    }
    auto_k <- TRUE
  } else {
    sel <- 1L
    auto_k <- FALSE
  }
  k_best <- scan$k[sel]
  fit <- best[[sel]]
  sil <- .emeans_silhouette(Z, fit$labels)

  ## ---- back to gene space: exact subtype means, no back-projection ----
  lab <- fit$labels
  sizes <- tabulate(lab, k_best)
  C_scaled <- matrix(0, k_best, ncol(prep$X))
  rs <- rowsum(prep$X, group = lab, reorder = TRUE)
  g <- as.integer(rownames(rs))
  C_scaled[g, ] <- rs / sizes[g]
  ## original scale: x = z * sd + mean
  C_expr <- sweep(C_scaled, 2, prep$scale, "*")
  C_expr <- sweep(C_expr, 2, prep$center, "+")
  rownames(C_scaled) <- rownames(C_expr) <- rownames(fit$centroids) <-
    paste0("C", seq_len(k_best))
  colnames(C_scaled) <- colnames(C_expr) <- prep$gene

  cluster <- factor(paste0("C", lab), levels = paste0("C", seq_len(k_best)))
  names(cluster) <- prep$sample

  if (verbose) {
    cat(sprintf("  selected k = %d%s   sizes: %s   silhouette %.3f\n",
                k_best, if (auto_k) " (auto)" else "",
                paste(sizes, collapse = " / "), sil$mean))
    cat(sprintf("  total %.1fs\n", proc.time()[["elapsed"]] - t_start))
  }

  structure(list(
    cluster        = cluster,
    k              = k_best,
    sizes          = sizes,
    mse            = fit$mse,
    silhouette     = sil$mean,
    silhouette_values = sil$values,
    silhouette_idx = sil$idx,
    k_scan         = scan,
    centroids      = C_scaled,
    centroids_expr = C_expr,
    centroids_pc   = fit$centroids,
    scores         = Z,
    pve            = prep$pve,
    fit            = fit,
    prep           = prep[c("gene", "sample", "center", "scale", "centered",
                            "scaled", "rotation", "npc", "gene_median")],
    expr_scaled    = prep$X,
    n_init         = as.integer(max(1L, n_init)),
    call           = match.call(),
    elapsed        = proc.time()[["elapsed"]] - t_start
  ), class = "emeans_bio")
}


#' @param x An `emeans_bio` fit.
#' @param ... Ignored.
#' @export
#' @rdname emeans_bio
print.emeans_bio <- function(x, ...) {
  cat(sprintf("<e-means subtype fit>  %d samples x %d genes  k = %d\n",
              length(x$cluster), ncol(x$expr_scaled), x$k))
  cat(sprintf("  sizes: %s\n",
              paste(sprintf("C%d=%d", seq_len(x$k), x$sizes), collapse = "  ")))
  cat(sprintf("  silhouette = %.3f   train MSE = %.4g   epochs = %d   restarts = %d   total %.1fs\n",
              x$silhouette, x$mse, x$fit$epochs, x$n_init,
              x$elapsed))
  cat("  use: plot(fit) | emeans_markers(fit) | predict(fit, new_expr)\n")
  invisible(x)
}


#' @param object An `emeans_bio` fit.
#' @param ... Ignored.
#' @export
#' @rdname emeans_bio
summary.emeans_bio <- function(object, ...) {
  cat(sprintf("e-means disease subtyping: %d samples, %d genes\n",
              length(object$cluster), ncol(object$expr_scaled)))
  cat(sprintf("selected k = %d, silhouette = %.3f\n\n",
              object$k, object$silhouette))
  print(object$k_scan, row.names = FALSE)
  cat("\ncluster sizes:\n"); print(table(object$cluster))
  invisible(object$k_scan)
}


.emeans_pal <- function(k) {
  pal <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
           "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD")
  rep(pal, length.out = max(k, 1))
}


#' Plot an e-means subtype fit
#'
#' Draws either the samples on the first two principal components, with the
#' subtype centroids marked by crosses, or the silhouette profile.
#'
#' @param x An `emeans_bio` fit.
#' @param type `"pca"` (default) or `"silhouette"`.
#' @param main Plot title; a sensible default is used when `NULL`.
#' @param ... Passed to [graphics::plot()].
#'
#' @return `x`, invisibly.
#'
#' @examples
#' \donttest{
#' fit <- emeans_bio(emeans_example, seed = 1, verbose = FALSE)
#' plot(fit)
#' plot(fit, type = "silhouette")
#' }
#'
#' @export
plot.emeans_bio <- function(x, type = c("pca", "silhouette"), main = NULL, ...) {
  type <- match.arg(type)
  cl <- x$cluster
  cols <- .emeans_pal(x$k)[as.integer(cl)]

  if (type == "pca") {
    ## when the fit skipped PCA, project on the fly; otherwise the first two
    ## genes would be drawn as if they were PC1/PC2
    if (!is.null(x$prep$rotation) && !is.na(x$prep$npc)) {
      Z  <- x$scores
      cp <- x$fit$centroids
      pv <- if (!is.null(x$pve)) 100 * x$pve[1:2] else c(NA, NA)
    } else {
      Xs <- x$expr_scaled
      mu <- colMeans(Xs)
      p2 <- prcomp(Xs, center = TRUE, scale. = FALSE, rank. = 2)
      Z  <- p2$x[, 1:2, drop = FALSE]
      cp <- sweep(x$fit$centroids, 2, mu, "-") %*% p2$rotation[, 1:2, drop = FALSE]
      tv <- sum(apply(Xs, 2, var))
      pv <- 100 * (p2$sdev[1:2]^2) / tv
    }
    if (is.null(Z) || ncol(Z) < 2) stop("not enough dimensions to plot", call. = FALSE)
    plot(Z[, 1], Z[, 2], col = cols, pch = 19, cex = 0.9,
         xlab = sprintf("PC1 (%.1f%%)", pv[1]),
         ylab = sprintf("PC2 (%.1f%%)", pv[2]),
         main = if (is.null(main)) sprintf("e-means subtypes (k = %d)", x$k) else main,
         ...)
    points(cp[, 1], cp[, 2], pch = 4, cex = 2.2, lwd = 2.5, col = "black")
    ## put the legend in the emptiest corner so it does not cover points
    xr <- range(Z[, 1]); yr <- range(Z[, 2])
    dx <- diff(xr); dy <- diff(yr)
    cand <- list(
      topright    = c(xr[2] - 0.35 * dx, xr[2], yr[2] - 0.35 * dy, yr[2]),
      topleft     = c(xr[1], xr[1] + 0.35 * dx, yr[2] - 0.35 * dy, yr[2]),
      bottomright = c(xr[2] - 0.35 * dx, xr[2], yr[1], yr[1] + 0.35 * dy),
      bottomleft  = c(xr[1], xr[1] + 0.35 * dx, yr[1], yr[1] + 0.35 * dy))
    cnt <- vapply(cand, function(b) {
      sum(Z[, 1] >= b[1] & Z[, 1] <= b[2] & Z[, 2] >= b[3] & Z[, 2] <= b[4])
    }, integer(1))
    pos <- names(which.min(cnt))[1]
    legend(pos, legend = levels(cl), col = .emeans_pal(x$k),
           pch = 19, bty = "n", title = "subtype")
  } else {
    s <- x$silhouette_values
    if (!length(s)) stop("silhouette values not stored", call. = FALSE)
    cls <- droplevels(cl[x$silhouette_idx])
    ord <- order(cls, -s)
    s <- s[ord]; cls <- cls[ord]
    cols2 <- .emeans_pal(x$k)[as.integer(cls)]
    barplot(s, col = cols2, border = NA, space = 0,
            ylim = c(min(0, min(s)), max(0.05, max(s)) * 1.12),
            main = if (is.null(main)) sprintf("silhouette (mean %.3f)", x$silhouette) else main,
            ylab = "silhouette", xlab = "samples sorted within subtype")
    abline(h = x$silhouette, lty = 2)
    abline(h = 0, col = "grey40")
  }
  invisible(x)
}
