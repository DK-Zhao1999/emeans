# emeans 0.1.0

First release, packaging the `e_means_bio.R` reference implementation as an
installable R package.

## Interface

* `emeans_bio()` takes a genes-by-samples matrix and returns subtype labels,
  expression centroids, silhouette diagnostics, the `k` scan and a predictor.
* `emeans()`, `emeans_law()` and `emeans_assign()` expose the clustering
  engine for matrices that are already samples by features.
* `emeans_markers()` reports marker genes per subtype.
* `print()`, `summary()`, `plot()` and `predict()` methods for `emeans_bio`
  objects, registered through `NAMESPACE` rather than by hand.

## Changes relative to the standalone script

* A centroid that receives no sample in a batch, with `alpha0 = 1`, used to
  be left at the origin: `N` drops to exactly zero and `S / N` is a vector of
  zeros, which passes an `is.finite()` check. Such centroids are now re-seeded
  from the batch. This only affects the `alpha0 = 1` case, which the scaling
  law selects for large `n`.
* `emeans_markers()` returns a zero-row data frame with the documented columns
  when every subtype has fewer than two samples, instead of `NULL`.
* Update history is preallocated instead of grown with `c()` at every update.
* The missing-value warning now counts the values that are actually imputed,
  not the ones in genes that are dropped for being empty.
* `emeans()` rejects `NA` input with an explicit message, since imputation is
  the job of `emeans_bio()`.
* Input validation for `k`, `epochs`, `alpha0`, `init` and `n` is explicit and
  returns informative errors.
