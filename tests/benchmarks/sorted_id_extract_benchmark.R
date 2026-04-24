# Optional benchmark: random vs sorted-unique (dense) vs sorted-unique (sparse)
# request orders for extract_sequences. Not run by testthat.
# Usage: Rscript tests/benchmarks/sorted_id_extract_benchmark.R
suppressPackageStartupMessages({
  library(fastqindexr)
})

set.seed(1L)
n_rec <- 10000L
k <- 2000L
tmp <- tempfile(fileext = ".fa.gz")
make_benchmark_fasta(tmp, n = n_rec, width = 120L)

cat("Building index ...\n")
print(system.time({
  idx <- create_index(tmp, type = "fasta")
}))

# Unsorted; duplicates possible (triggers per-file sort + unique).
ids_random <- sample.int(n_rec, k, replace = TRUE)
# Sorted unique: dense local ID ranges (skips per-file sort in C++ when global
# order is sorted and unique).
ids_sorted_dense <- sort(sample.int(n_rec, k, replace = FALSE))
# Evenly spaced: few records per region, many regions.
ids_sorted_sparse <- as.integer(seq(1L, n_rec, length.out = k))

cat("\nRandom / unsorted IDs (with replacement)\n")
print(system.time({
  r1 <- extract_sequences(idx, ids_random)
}))

cat("\nSorted unique (dense line ranges expected)\n")
print(system.time({
  r2 <- extract_sequences(idx, ids_sorted_dense)
}))

cat("\nSorted unique (evenly spaced, sparse in line index)\n")
print(system.time({
  r3 <- extract_sequences(idx, ids_sorted_sparse)
}))

cat("\nRow counts:", nrow(r1), nrow(r2), nrow(r3), "\n")
unlink(tmp)
