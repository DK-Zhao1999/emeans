## -----------------------------------------------------------------------------
## Preprocessing: genes x samples -> samples x genes, NA imputation,
## zero-variance filter, gene-wise z-score, PCA
## -----------------------------------------------------------------------------
.emeans_prep <- function(x, center = TRUE, scale = TRUE, pca = 50) {
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("`x` must be a numeric matrix (genes x samples)", call. = FALSE)
  }
  if (nrow(x) < 2 || ncol(x) < 2) {
    stop("`x` needs at least 2 genes and 2 samples", call. = FALSE)
  }
  ## Only warn about orientation when the matrix is wide enough that it really
  ## looks transposed: a small gene panel measured on thousands of samples is
  ## a legitimate input.
  if (nrow(x) < ncol(x) && ncol(x) >= 2000L) {
    warning("ncol(x) is much larger than nrow(x), but this function expects ",
            "genes in rows and samples in columns; check the orientation",
            call. = FALSE)
  }

  gene <- rownames(x)
  samp <- colnames(x)
  X <- t(x)                                   # samples x genes
  storage.mode(X) <- "double"
  if (is.null(gene)) gene <- sprintf("G%05d", seq_len(ncol(X)))
  if (is.null(samp)) samp <- sprintf("S%04d", seq_len(nrow(X)))

  ## NA: impute the gene median; genes that are missing everywhere are dropped
  n_na <- 0L
  if (anyNA(X)) {
    allna <- which(apply(is.na(X), 2, all))
    if (length(allna)) {
      warning(sprintf("%d genes are entirely missing and were dropped",
                      length(allna)), call. = FALSE)
      X <- X[, -allna, drop = FALSE]
      gene <- gene[-allna]
      if (ncol(X) < 2) {
        stop("fewer than 2 genes with data left; check the input", call. = FALSE)
      }
    }
    n_na <- sum(is.na(X))
    if (n_na > 0L) {
      med <- apply(X, 2, median, na.rm = TRUE)
      for (j in seq_len(ncol(X))) {
        bad <- is.na(X[, j])
        if (any(bad)) X[bad, j] <- med[j]
      }
      warning(sprintf("%d missing values imputed with gene medians", n_na),
              call. = FALSE)
    }
  }
  ## Inf/NaN are not handled silently: a gene slipping away unnoticed is worse
  ## than a hard stop
  if (any(!is.finite(X))) {
    stop(sprintf(paste0("`x` contains %d non-finite values (Inf/-Inf); ",
                        "please remove or fix them first"),
                 sum(!is.finite(X))), call. = FALSE)
  }
  gene_median <- apply(X, 2, median)

  mu <- colMeans(X)
  sdv <- apply(X, 2, sd)
  keep <- is.finite(sdv) & sdv > 0
  if (!all(keep)) {
    warning(sprintf("%d zero-variance genes dropped", sum(!keep)), call. = FALSE)
    X <- X[, keep, drop = FALSE]; gene <- gene[keep]
    mu <- mu[keep]; sdv <- sdv[keep]; gene_median <- gene_median[keep]
  }
  if (ncol(X) < 2) {
    stop("fewer than 2 non-constant genes left; check the input", call. = FALSE)
  }

  Xs <- X
  if (center) Xs <- sweep(Xs, 2, mu, "-")
  if (scale)  Xs <- sweep(Xs, 2, sdv, "/")

  out <- list(
    X = Xs, gene = gene, sample = samp,
    center = if (center) mu else rep(0, ncol(Xs)),
    scale  = if (scale) sdv else rep(1, ncol(Xs)),
    centered = center, scaled = scale,
    gene_median = gene_median, n_na = n_na,
    npc = NA_integer_, rotation = NULL, scores = Xs
  )

  ## PCA: with far fewer samples than genes the Gram matrix plus eigen() is an
  ## order of magnitude faster than prcomp; both give the same scores and
  ## rotation, up to the sign of each component
  if (!is.null(pca) && is.finite(pca) && pca > 0) {
    npc <- as.integer(min(pca, nrow(Xs) - 1L, ncol(Xs)))
    tot_ss <- sum(Xs * Xs)
    if (npc >= 2) {
      if (nrow(Xs) <= ncol(Xs) && nrow(Xs) <= 1500L) {
        eg <- eigen(tcrossprod(Xs), symmetric = TRUE)
        lam <- pmax(eg$values, 0)
        rk <- sum(lam > max(lam) * 1e-10)
        npc <- if (rk >= 2L) min(npc, rk) else 0L
        if (npc >= 2L) {
          U <- eg$vectors[, seq_len(npc), drop = FALSE]
          sv <- sqrt(lam[seq_len(npc)])
          scores <- sweep(U, 2, sv, "*")
          rotation <- crossprod(Xs, sweep(U, 2, 1 / sv, "*"))
          pve <- lam[seq_len(npc)] / max(tot_ss, .Machine$double.eps)
        }
      } else {
        pc <- prcomp(Xs, center = FALSE, scale. = FALSE, rank. = npc)
        scores <- pc$x[, seq_len(npc), drop = FALSE]
        rotation <- pc$rotation[, seq_len(npc), drop = FALSE]
        pve <- (pc$sdev^2 * (nrow(Xs) - 1)) / max(tot_ss, .Machine$double.eps)
        pve <- pve[seq_len(npc)]
      }
      if (npc >= 2L) {
        colnames(scores) <- paste0("PC", seq_len(npc))
        colnames(rotation) <- paste0("PC", seq_len(npc))
        rownames(rotation) <- colnames(Xs)
        out$rotation <- rotation
        out$scores <- scores
        out$npc <- npc
        out$pve <- pve
      }
    }
  }
  out
}


## -----------------------------------------------------------------------------
## Mean silhouette width, without the cluster package.
## Above `max_n` samples the width is estimated from an evenly spaced subsample
## so that the n x n distance matrix cannot exhaust memory.
## -----------------------------------------------------------------------------
.emeans_silhouette <- function(Z, cl, max_n = 2000L) {
  n <- nrow(Z)
  if (n < 3) return(list(mean = NA_real_, values = numeric(0), idx = integer(0)))
  idx <- if (n > max_n) unique(round(seq(1, n, length.out = max_n))) else seq_len(n)
  cl <- cl[idx]
  lev <- sort(unique(cl))
  K <- length(lev)
  if (K < 2 || K >= length(cl)) {
    return(list(mean = NA_real_, values = numeric(0), idx = idx))
  }
  d <- as.matrix(dist(Z[idx, , drop = FALSE]))
  s <- numeric(length(cl))
  for (i in seq_along(cl)) {
    own <- which(cl == cl[i])
    if (length(own) <= 1L) { s[i] <- 0; next }
    a <- mean(d[i, own[own != i]])
    b <- min(vapply(lev[lev != cl[i]], function(g) mean(d[i, cl == g]), numeric(1)))
    s[i] <- if (max(a, b) > 0) (b - a) / max(a, b) else 0
  }
  list(mean = mean(s), values = s, idx = idx)
}
