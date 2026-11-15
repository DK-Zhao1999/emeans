test_that("emeans_law returns a usable schedule", {
  for (n in c(50, 500, 5000)) {
    for (k in c(2, 5)) {
      l <- emeans_law(n, k)
      expect_true(l$batch >= 16 || l$batch == n)
      expect_lte(l$batch, n)
      expect_equal(log2(l$batch), round(log2(l$batch)))   # power of two
      expect_gt(l$alpha0, 0)
      expect_lte(l$alpha0, 1)
    }
  }
  ## small n: the batch is clamped to n
  expect_equal(emeans_law(10, 2)$batch, 10)
})

test_that("emeans separates well-separated groups", {
  set.seed(1)
  x <- rbind(matrix(rnorm(100 * 2, 0), 100, 2),
             matrix(rnorm(100 * 2, 6), 100, 2),
             matrix(rnorm(100 * 2, -6), 100, 2))
  truth <- rep(1:3, each = 100)
  fit <- emeans(x, k = 3, seed = 1)
  expect_s3_class(fit, "emeans")
  expect_equal(sort(tabulate(fit$labels, 3)), c(100, 100, 100))
  ## every recovered cluster is one true group
  expect_equal(length(unique(paste(truth, fit$labels))), 3)
  expect_length(fit$history$update, fit$updates)
})

test_that("emeans validates its input", {
  x <- matrix(rnorm(40), 20, 2)
  expect_error(emeans(x, k = 0), "positive integer")
  expect_error(emeans(x, k = 2, epochs = 0), "epochs")
  expect_error(emeans(`colnames<-`(x, NULL) + NA_real_, k = 2), "NA")
  expect_warning(emeans(x, k = 50, seed = 1), "k > n")
})

test_that("emeans handles an empty cluster without collapsing a centroid", {
  set.seed(3)
  x <- matrix(rnorm(120 * 3, mean = 5), 120, 3)
  fit <- emeans(x, k = 8, epochs = 3, batch = 4, alpha0 = 1,
                revival = FALSE, seed = 1)
  ## a centroid that receives no sample must be re-seeded from the data, not
  ## left at the origin
  norms <- sqrt(rowSums(fit$centroids^2))
  expect_true(all(norms > 0))
})

test_that("emeans_bio recovers the simulated subtypes", {
  fit <- emeans_bio(emeans_example, k = 3, seed = 1, verbose = FALSE)
  expect_s3_class(fit, "emeans_bio")
  expect_equal(fit$k, 3)
  expect_equal(fit$sizes, c(40, 40, 40))
  truth <- rep(c("S1", "S2", "S3"), each = 40)
  ## cross-tabulation is a permutation matrix: each truth group lands in one fit
  tab <- table(truth, fit$cluster)
  expect_true(all(rowSums(tab > 0) == 1))
  expect_true(all(colSums(tab > 0) == 1))
  expect_gt(fit$silhouette, 0.2)
})

test_that("emeans_bio selects k automatically", {
  fit <- emeans_bio(emeans_example, k_max = 4, n_init = 2, epochs = 3,
                    seed = 1, verbose = FALSE)
  expect_true(fit$k %in% 2:4)
  expect_equal(nrow(fit$k_scan), 3L)
  expect_true(all(c("k", "silhouette", "mse", "seconds") %in% names(fit$k_scan)))
})

test_that("emeans_bio handles a single fixed k and no PCA", {
  fit <- emeans_bio(emeans_example, k = 3, pca = NULL, epochs = 2,
                    n_init = 1, seed = 2, verbose = FALSE)
  expect_true(is.na(fit$prep$npc))
  expect_null(fit$prep$rotation)
  expect_equal(ncol(fit$scores), nrow(emeans_example))
  expect_equal(nrow(fit$centroids_expr), 3L)
})

test_that("centroids are exact subtype means on both scales", {
  fit <- emeans_bio(emeans_example, k = 3, seed = 1, verbose = FALSE)
  ## in the scaled space the centroid is the exact mean of its samples
  for (j in seq_len(fit$k)) {
    direct <- colMeans(fit$expr_scaled[fit$cluster == levels(fit$cluster)[j], ,
                                       drop = FALSE])
    expect_equal(unname(fit$centroids[j, ]), unname(direct), tolerance = 1e-10)
  }
  ## and centroids_expr is that centroid carried back to the input scale
  back <- sweep(fit$centroids, 2, fit$prep$scale, "*")
  back <- sweep(back, 2, fit$prep$center, "+")
  expect_equal(unname(fit$centroids_expr), unname(back), tolerance = 1e-12)
})

test_that("input problems are reported clearly", {
  expect_error(emeans_bio(matrix(1, 3, 3)), "zero-variance|non-constant")
  bad <- emeans_example
  bad[1, 1] <- Inf
  expect_error(emeans_bio(bad, k = 3, verbose = FALSE), "non-finite")
  nas <- emeans_example
  nas[1:3, 1:3] <- NA
  expect_warning(emeans_bio(nas, k = 3, n_init = 1, epochs = 2, seed = 1,
                            verbose = FALSE), "imputed")
})

test_that("markers are ranked and typed", {
  fit <- emeans_bio(emeans_example, k = 3, seed = 1, verbose = FALSE)
  m <- emeans_markers(fit, n = 5)
  expect_s3_class(m, "data.frame")
  expect_equal(nrow(m), 15L)
  expect_true(all(c("cluster", "gene", "mean_in", "mean_out", "diff",
                    "t", "p", "padj") %in% names(m)))
  expect_true(all(m$p >= 0 & m$p <= 1))
  expect_false(is.unsorted(m$p[m$cluster == m$cluster[1]]))
})

test_that("markers return an empty frame when every subtype is too small", {
  small <- emeans_bio(emeans_example[, 1:3, drop = FALSE], k = 2, epochs = 1,
                      n_init = 1, seed = 3, verbose = FALSE)
  m <- emeans_markers(small, n = 2)
  expect_s3_class(m, "data.frame")
  expect_equal(nrow(m), 0L)
  expect_true("padj" %in% names(m))
})

test_that("predict reproduces the training labels", {
  fit <- emeans_bio(emeans_example, k = 3, seed = 1, verbose = FALSE)
  pred <- predict(fit, emeans_example)
  expect_s3_class(pred, "factor")
  expect_equal(names(pred), colnames(emeans_example))
  expect_equal(as.character(pred), as.character(fit$cluster))
})

test_that("predict works on a single sample and on reordered genes", {
  fit <- emeans_bio(emeans_example, k = 3, seed = 1, verbose = FALSE)
  one <- emeans_example[, 1, drop = FALSE]
  expect_length(predict(fit, one), 1L)

  shuffled <- emeans_example[sample(nrow(emeans_example)), , drop = FALSE]
  expect_equal(as.character(predict(fit, shuffled)),
               as.character(predict(fit, emeans_example)))

  sc <- predict(fit, one, type = "score")
  expect_true(all(c("cluster", "dist", "score") %in% names(sc)))
})

test_that("predict refuses a matrix with missing genes", {
  fit <- emeans_bio(emeans_example, k = 3, seed = 1, verbose = FALSE)
  expect_error(predict(fit, emeans_example[-1, , drop = FALSE]), "missing")
})

test_that("print, summary and plot run", {
  fit <- emeans_bio(emeans_example, k = 3, seed = 1, verbose = FALSE)
  expect_output(print(fit), "e-means subtype fit")
  expect_output(summary(fit), "cluster sizes")
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_invisible(plot(fit))
  expect_invisible(plot(fit, type = "silhouette"))
})
