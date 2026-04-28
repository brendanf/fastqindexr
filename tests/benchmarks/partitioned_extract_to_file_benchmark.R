#!/usr/bin/env Rscript

pkgload::load_all(".", quiet = TRUE)

tmp <- tempfile(fileext = ".fa.gz")
on.exit(unlink(tmp), add = TRUE)
make_benchmark_fasta(tmp, n = 250000L, width = 60L)
idx <- create_index(tmp, type = "fasta")

ids <- seq_len(200000L)
cont <- partition_seq_idx(ids, n_parts = 8L, strategy = "contiguous")
rr <- partition_seq_idx(ids, n_parts = 8L, strategy = "round_robin")

bench_partition <- function(label, partitions) {
  outs <- replicate(
    length(partitions),
    tempfile(fileext = ".fa"),
    simplify = TRUE
  )
  on.exit(unlink(outs), add = TRUE)
  elapsed <- system.time(
    extract_sequences_to_file(idx, partitions, outfile = outs, append = FALSE)
  )[["elapsed"]]
  cat(sprintf("%s elapsed: %.3fs\n", label, elapsed))
}

cat("Partitioned strict-increasing extraction benchmark\n")
bench_partition("contiguous-8", cont)
bench_partition("round_robin-8", rr)
