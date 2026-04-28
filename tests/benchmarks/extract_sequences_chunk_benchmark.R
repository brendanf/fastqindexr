#!/usr/bin/env Rscript
# Benchmark `extract_sequences()` after the C++ backend processes unique record
# IDs in fixed-size chunks (see `kChunkUniqueIds` in `cpp_extract_sequences`).
# This is separate from `extract_sequences_dnastringset(..., chunk_chars = )`.
# That argument is an R-visible tuning knob: it must be a finite positive scalar
# (`NA` and `Inf` are errors). Run from package root:
# Rscript tests/benchmarks/extract_sequences_chunk_benchmark.R

pkgload::load_all(".", quiet = TRUE)

bench_case <- function(label, idx, ids, return = "data.frame") {
  gc(reset = TRUE)
  t <- system.time({
    out <- extract_sequences(idx, ids, return = return)
  })
  cat(sprintf(
    "%s return=%s elapsed=%.3fs n=%d\n",
    label,
    return,
    t[["elapsed"]],
    length(ids)
  ))
  invisible(out)
}

set.seed(1L)
tmp <- tempfile(fileext = ".fa.gz")
on.exit(unlink(tmp), add = TRUE)
make_benchmark_fasta(tmp, n = 500000L, width = 120L)
idx <- create_index(tmp, type = "fasta")

cat("extract_sequences chunked-unique backend benchmark\n")

ordered_unique <- seq_len(300000L)
bench_case("ordered_unique", idx, ordered_unique, return = "data.frame")
bench_case("ordered_unique", idx, ordered_unique, return = "seq")

set.seed(2L)
random_dup <- sample.int(500000L, 300000L, replace = TRUE)
bench_case("random_with_replacement", idx, random_dup, return = "data.frame")
bench_case("random_with_replacement", idx, random_dup, return = "seq")

set.seed(3L)
perm_unique <- sample.int(500000L, 300000L, replace = FALSE)
bench_case("permuted_unique", idx, perm_unique, return = "data.frame")
bench_case("permuted_unique", idx, perm_unique, return = "seq")
