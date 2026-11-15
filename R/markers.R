#' Marker genes for each subtype
#'
#' For every subtype, compares its samples against all remaining samples gene
#' by gene with a Welch t-test on the scaled expression and returns the `n`
#' most significant genes.
#'
#' The test is a screening heuristic, not a differential-expression pipeline:
#' samples are not independent of the clustering, so the p-values are
#' optimistic. Use them to rank genes, and confirm candidates on an
#' independent cohort.
#'
#' @param fit An `emeans_bio` fit.
#' @param n Number of genes reported per subtype, default 20.
#'
#' @return A data frame with one row per gene and subtype, columns `cluster`,
#'   `gene`, `mean_in`, `mean_out`, `diff` (difference of scaled means), `t`,
#'   `p` and `padj` (Benjamini-Hochberg across genes within the subtype), sorted
#'   by `p`. Subtypes with fewer than two samples in either group contribute no
#'   rows; if none contribute, an empty data frame with the same columns is
#'   returned.
#'
#' @examples
#' fit <- emeans_bio(emeans_example, k = 3, seed = 1, verbose = FALSE)
#' markers <- emeans_markers(fit, n = 5)
#' head(markers)
#'
#' @seealso [emeans_bio()]
#'
#' @export
emeans_markers <- function(fit, n = 20L) {
  if (!inherits(fit, "emeans_bio")) {
    stop("`fit` must come from emeans_bio()", call. = FALSE)
  }
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("`n` must be a positive integer", call. = FALSE)
  }
  X <- fit$expr_scaled                 # samples x genes, already scaled
  cl <- fit$cluster
  lev <- levels(cl)
  out <- vector("list", length(lev))

  for (j in seq_along(lev)) {
    g1 <- X[cl == lev[j], , drop = FALSE]
    g2 <- X[cl != lev[j], , drop = FALSE]
    n1 <- nrow(g1); n2 <- nrow(g2)
    if (n1 < 2L || n2 < 2L) next
    m1 <- colMeans(g1); m2 <- colMeans(g2)
    v1 <- apply(g1, 2, var); v2 <- apply(g2, 2, var)
    se <- sqrt(v1 / n1 + v2 / n2)
    se[!is.finite(se) | se == 0] <- NA_real_
    tt <- (m1 - m2) / se
    df <- (v1 / n1 + v2 / n2)^2 /
      ((v1 / n1)^2 / (n1 - 1) + (v2 / n2)^2 / (n2 - 1))
    p <- 2 * pt(-abs(tt), df)
    p[!is.finite(p)] <- 1

    dfj <- data.frame(cluster = lev[j], gene = fit$prep$gene,
                      mean_in = m1, mean_out = m2, diff = m1 - m2,
                      t = tt, p = p, stringsAsFactors = FALSE)
    dfj$padj <- p.adjust(dfj$p, method = "BH")
    dfj <- dfj[order(dfj$p), , drop = FALSE]
    out[[j]] <- head(dfj, n)
  }
  res <- do.call(rbind, out)
  if (is.null(res)) {
    res <- data.frame(cluster = character(0), gene = character(0),
                      mean_in = numeric(0), mean_out = numeric(0),
                      diff = numeric(0), t = numeric(0), p = numeric(0),
                      padj = numeric(0), stringsAsFactors = FALSE)
  }
  rownames(res) <- NULL
  res
}
