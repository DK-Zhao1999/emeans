## Build data/emeans_example.rda: a simulated 3-subtype bulk expression matrix.
set.seed(20260101)

n_gene <- 300L
per    <- 40L
n_samp <- per * 3L
subtype <- rep(c("S1", "S2", "S3"), each = per)

genes   <- sprintf("G%04d", seq_len(n_gene))
samples <- sprintf("S%03d", seq_len(n_samp))

## baseline expression on a log2 scale; genes differ a lot in level
base <- rnorm(n_gene, mean = 8, sd = 1.5)

## three disjoint blocks of 40 marker genes, each raised by 2 in its subtype
shift <- matrix(0, n_gene, n_samp)
blocks <- list(1:40, 41:80, 81:120)
for (s in seq_len(3)) {
  shift[blocks[[s]], subtype == paste0("S", s)] <- 2
}

## shared latent factors give genes some correlation, as in real data, but are
## kept small enough that the subtypes stay separable
n_fac <- 6
fac <- matrix(rnorm(n_fac * n_samp, sd = 0.5), n_fac, n_samp)
loading <- matrix(0, n_gene, n_fac)
for (f in seq_len(n_fac)) {
  block <- ((f - 1) * 50 + 1):(f * 50)
  loading[block, f] <- runif(length(block), 0.2, 0.6)
}

expr <- base + shift + loading %*% fac +
  matrix(rnorm(n_gene * n_samp, sd = 0.4), n_gene, n_samp)
rownames(expr) <- genes
colnames(expr) <- samples

## a few missing values so the imputation path is exercised by the examples
na_idx <- cbind(sample(n_gene, 8), sample(n_samp, 8))
expr[na_idx] <- NA

emeans_example <- expr

dir.create("data", showWarnings = FALSE)
save(emeans_example, file = "data/emeans_example.rda", compress = "xz", version = 2)

cat("dim:", paste(dim(emeans_example), collapse = " x "),
    " NA:", sum(is.na(emeans_example)),
    " range:", paste(round(range(emeans_example, na.rm = TRUE), 1), collapse = " .. "), "\n")
cat("file:", round(file.size("data/emeans_example.rda") / 1024), "KB\n")

## sanity: does the documented pipeline recover the three blocks?
e <- new.env(parent = globalenv())
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) sys.source(f, envir = e)
fit <- e$emeans_bio(emeans_example, seed = 1, verbose = FALSE)
cat("recovered k =", fit$k, " sizes =", paste(fit$sizes, collapse = "/"),
    " silhouette =", round(fit$silhouette, 3), "\n")
cat("contingency (truth x fit):\n")
print(table(subtype, fit$cluster))
