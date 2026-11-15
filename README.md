# emeans

<!-- badges: start -->
[![R-CMD-check](https://github.com/DK-Zhao1999/emeans/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DK-Zhao1999/emeans/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
<!-- badges: end -->

**Unsupervised subtype discovery for bulk expression data, in one call.**

`emeans` clusters samples with [e-means](https://e-means.github.io), a k-means
variant whose centroids follow an exponentially decaying average, with a cosine
annealed step size and periodic revival of dead centroids. You give it a
genes-by-samples matrix; it gives back subtype labels, expression centroids,
marker genes and a predictor for new samples.

No Bioconductor dependency, no compiled code, no preprocessing pipeline to
assemble: median imputation, zero-variance filtering, gene-wise z-scores and
PCA are handled inside.

## Installation

```r
# install.packages("remotes")
remotes::install_github("DK-Zhao1999/emeans")
```

Requires R >= 3.5.0. The only imports are `stats`, `graphics` and `utils`.

## Quick start

```r
library(emeans)

expr <- emeans_example          # 300 genes x 120 samples, simulated, 3 subtypes
fit <- emeans_bio(expr)
#>   scanning k = 2-3-4-5-6 (silhouette), 10 inits each
#>   k=2  silhouette 0.191  train MSE 199.5  (0.1s)
#>   k=3  silhouette 0.280  train MSE 149.7  (0.1s)
#>   k=4  silhouette 0.206  train MSE 145.0  (0.1s)
#>   ...
#>   selected k = 3 (auto)   sizes: 40 / 40 / 40   silhouette 0.280

fit
#> <e-means subtype fit>  120 samples x 300 genes  k = 3
#>   sizes: C1=40  C2=40  C3=40
#>   silhouette = 0.280   train MSE = 149.7   epochs = 5   restarts = 10   total 0.6s
```

That is the whole workflow. From there:

```r
summary(fit)                    # the k scan: silhouette / error / time per k
plot(fit)                       # PC1-PC2, centroids marked with crosses
plot(fit, type = "silhouette")  # silhouette profile

head(emeans_markers(fit), 5)    # marker genes per subtype
predict(fit, new_expr)          # assign new samples to the learned subtypes
```

## What you get back

| Element | What it is |
| --- | --- |
| `cluster` | factor of subtype labels, named by sample |
| `sizes` | samples per subtype |
| `centroids_expr` | subtype by gene matrix of means, **on the input scale** |
| `centroids` | the same centroids in the scaled space used for clustering |
| `silhouette`, `k_scan` | diagnostics, including the full `k` scan |
| `scores`, `pve` | sample by PC coordinates and variance explained |
| `fit` | the underlying `emeans` engine object |

## How it works

1. Transpose to samples by genes, impute missing values with the gene median,
   drop genes that are empty everywhere and genes with zero variance.
2. z-score each gene.
3. Keep the first `pca` (default 50) principal components.
4. Run e-means: batch-wise updates of `S` and `N` with a cosine-annealed step
   size, splitting of under-used centroids, and a final pass that recomputes
   each centroid as the exact mean of the samples it owns.
5. Choose `k` by mean silhouette width over `2:k_max` (default 6), unless you
   pass `k` yourself. Each candidate `k` is run `n_init` (default 10) times and
   the lowest-error run is kept.
6. Recompute the centroids in the gene space, as subtype means.

The batch size and initial step size come from the scaling law in the paper
(`emeans_law()`),
so the number of updates stays roughly constant as the cohort grows.

## Choosing k yourself

```r
fit <- emeans_bio(expr, k = 4)
```

Automatic selection is a convenience, not a verdict: silhouette differences
between neighbouring `k` can be small on real cohorts. Look at `fit$k_scan`,
`plot(fit)` and the biological coherence of `emeans_markers(fit)` before
committing to a number.

## Notes and limits

- Clustering happens in PC space, so `silhouette` is measured there too and is
  not comparable with a silhouette computed on raw expression distances.
- Marker p-values come from a Welch t-test of each subtype against the rest.
  Samples are not independent of the clustering, so the values are optimistic;
  use them to rank genes and validate candidates on independent data.
- Above 2000 samples the silhouette is estimated from an evenly spaced
  subsample, which keeps the distance matrix from exhausting memory.
- Built for bulk data: dozens to thousands of samples, thousands of genes. For
  millions of cells, use the reference faiss implementation or pseudobulk first.
- `pca = NULL` clusters on all genes; `finalize = "slow-ema"` uses a
  longer-running average and tends to generalize better to new samples.

## Citation

The algorithm is described in:

> Bessa, Mazaré, Szilvasy, Touvron, Picard, Douze, Jégou (2026). *e-means:
> outperforming k-means (in one-tenth the epochs)*. <https://e-means.github.io>

If you use the package, please cite that work, and the package itself:

```r
citation("emeans")
```

## Attribution

This package is an R implementation of the e-means algorithm. It follows the
reference implementation published by **Meta Platforms, Inc. and affiliates**
at <https://github.com/e-means/e-means.github.io>, which is released under the
MIT license. That copyright notice is retained in
[LICENSE.md](LICENSE.md); the package is distributed under the same MIT terms.

## License

MIT. See [LICENSE.md](LICENSE.md). Copyright (c) 2026 Dingkang Zhao and
Copyright (c) Meta Platforms, Inc. and affiliates.

