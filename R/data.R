#' Simulated bulk expression matrix with three subtypes
#'
#' A genes-by-samples matrix simulated to mimic a bulk transcriptomic cohort
#' with three molecular subtypes, small enough to run in a second and used by
#' the examples and tests.
#'
#' Simulated as \eqn{8 + \text{subtype shift} + \text{latent factors} +
#' \text{noise}} on a log2 scale: three disjoint blocks of 40 genes each carry a
#' shift of +2 in one subtype, six latent factors shared by blocks of 50 genes
#' introduce gene-gene correlation, and the noise standard deviation is 0.4.
#' There are no missing values; the imputation path is exercised by the tests.
#'
#' @format A numeric matrix with 300 rows (genes) and 120 columns (samples),
#'   samples ordered by subtype in three consecutive blocks of 40. Row names are
#'   gene identifiers (`G0001`-`G0300`), column names are sample identifiers
#'   (`S001`-`S120`). The generating script is in `data-raw/make_example.R`.
#'
#' @source Simulated, no external data.
#'
#' @examples
#' dim(emeans_example)
#' table(rep(c("S1", "S2", "S3"), each = 40))
"emeans_example"
