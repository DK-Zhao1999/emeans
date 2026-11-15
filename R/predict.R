#' Assign new samples to the discovered subtypes
#'
#' Applies the preprocessing stored in the fit (gene centring and scaling, the
#' PCA rotation) to a new expression matrix and assigns each sample to the
#' nearest subtype centroid.
#'
#' `newdata` must contain every gene the model was trained on, with the same
#' names. Extra genes are ignored; missing ones are an error, because silently
#' shifting a sample in PC space would produce a confident but meaningless
#' label. Missing values are imputed with the training medians.
#'
#' @param object An `emeans_bio` fit.
#' @param newdata Numeric matrix or data frame, genes in rows and samples in
#'   columns, in the same orientation as the training data.
#' @param type `"cluster"` (default) returns the labels, `"score"` returns a
#'   list with the labels, the distance to the assigned centroid and the PC
#'   coordinates.
#' @param ... Ignored.
#'
#' @return A factor of subtype labels named after the samples, or, for
#'   `type = "score"`, a list with `cluster`, `dist` and `score`.
#'
#' @examples
#' fit <- emeans_bio(emeans_example, k = 3, seed = 1, verbose = FALSE)
#' new_samples <- emeans_example[, 1:4, drop = FALSE]
#' predict(fit, new_samples)
#'
#' @seealso [emeans_bio()]
#'
#' @export
predict.emeans_bio <- function(object, newdata, type = c("cluster", "score"), ...) {
  type <- match.arg(type)
  Xn <- as.matrix(newdata)
  if (!is.numeric(Xn)) {
    stop("`newdata` must be numeric (genes x samples)", call. = FALSE)
  }
  Xn <- t(Xn)
  if (!is.matrix(Xn)) Xn <- matrix(Xn, nrow = 1L)

  ## put the genes in the training order
  gn <- rownames(newdata)
  if (!is.null(gn)) {
    idx <- match(object$prep$gene, gn)
    if (anyNA(idx)) {
      stop(sprintf("newdata is missing %d genes used in training (e.g. %s)",
                   sum(is.na(idx)),
                   paste(utils::head(gn[is.na(idx)]), collapse = ", ")),
           call. = FALSE)
    }
    Xn <- Xn[, idx, drop = FALSE]
  } else {
    if (ncol(Xn) != length(object$prep$gene)) {
      stop("newdata has no rownames and its gene count differs from training",
           call. = FALSE)
    }
    warning("newdata has no rownames: assuming genes are in training order",
            call. = FALSE)
  }
  ## NA: impute with the training medians
  if (anyNA(Xn)) {
    for (j in seq_len(ncol(Xn))) {
      bad <- is.na(Xn[, j])
      if (any(bad)) Xn[bad, j] <- object$prep$gene_median[j]
    }
  }
  Xn <- sweep(Xn, 2, object$prep$center, "-")
  if (isTRUE(object$prep$scaled)) Xn <- sweep(Xn, 2, object$prep$scale, "/")

  Z <- if (!is.null(object$prep$rotation)) Xn %*% object$prep$rotation else Xn
  asg <- emeans_assign(Z, object$centroids_pc)

  out <- factor(levels(object$cluster)[asg$labels], levels = levels(object$cluster))
  names(out) <- if (!is.null(colnames(newdata))) colnames(newdata) else
    sprintf("S%04d", seq_len(nrow(Z)))
  if (type == "score") {
    return(list(cluster = out, dist = asg$dist, score = Z))
  }
  out
}
