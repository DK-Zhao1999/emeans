#' emeans: fast unsupervised subtyping of expression data
#'
#' `emeans` implements the e-means clustering algorithm and wraps it in a
#' single call that turns a genes-by-samples expression matrix into per-sample
#' subtype labels.
#'
#' e-means is a k-means variant. Instead of moving each centroid to the mean of
#' the points it currently owns, it keeps a running average
#' (\eqn{S_j / N_j}) that is updated with an exponentially decaying step size,
#' annealed with a cosine schedule. A revival step splits under-used centroids
#' so the search does not stall in a local optimum, and a final pass recomputes
#' the centroids from the last epoch's assignment.
#'
#' @section Typical use:
#'
#' ```r
#' fit <- emeans_bio(expr)      # expr: genes x samples
#' fit                          # selected k and cluster sizes
#' plot(fit)                    # PC1-PC2 with centroids
#' emeans_markers(fit)          # marker genes per subtype
#' predict(fit, new_expr)       # assign new samples
#' ```
#'
#' @section Where to start:
#'
#' * [emeans_bio()] is the interface for expression data.
#' * [emeans()] is the clustering engine, for matrices that are already
#'   samples by features.
#' * [emeans_markers()] and [predict.emeans_bio()] work on a fitted object.
#' * [emeans_law()] documents the batch/step-size schedule.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom graphics abline barplot legend points
#' @importFrom stats dist median p.adjust prcomp pt rnorm sd var
#' @importFrom utils head
NULL
