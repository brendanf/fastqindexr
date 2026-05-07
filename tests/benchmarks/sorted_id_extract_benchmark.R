# Optional benchmark: random vs sorted-unique (dense) vs sorted-unique (sparse)
# request orders for extract_sequences. Not run by testthat.
# Usage: Rscript tests/benchmarks/sorted_id_extract_benchmark.R
pkgload::load_all(".", quiet = TRUE)

source(file.path("tests", "benchmarks", "benchmark_sink.R"))
bench_sink_start("sorted_id_extract_benchmark")
on.exit(bench_sink_stop(), add = TRUE)

set.seed(1L)
n_rec <- 10000L
k <- 2000L

run_one <- function(fmt) {
  ext <- if (fmt == "fasta") ".fa.gz" else ".fq.gz"
  tmp <- tempfile(fileext = ext)
  on.exit(unlink(tmp), add = TRUE)
  if (fmt == "fasta") {
    make_benchmark_fasta(tmp, n = n_rec, width = 120L)
  } else {
    make_benchmark_fastq(tmp, n = n_rec, width = 120L)
  }

  cat(sprintf("\n[benchmark]\n%s: build index\n", toupper(fmt)))
  t_idx <- system.time({
    idx <- create_index(tmp, type = fmt)
  })
  cat(sprintf("case: %s/index\n", fmt))
  cat(sprintf("elapsed_s: %.3f\n", t_idx[["elapsed"]]))

  # Unsorted; duplicates possible (triggers per-file sort + unique).
  ids_random <- sample.int(n_rec, k, replace = TRUE)
  # Sorted unique: dense local ID ranges.
  ids_sorted_dense <- sort(sample.int(n_rec, k, replace = FALSE))
  # Evenly spaced: sparse in local ID space.
  ids_sorted_sparse <- as.integer(seq(1L, n_rec, length.out = k))

  t1 <- system.time({
    r1 <- extract_sequences(idx, ids_random)
  })
  cat(sprintf("case: %s/random_unsorted_dup\n", fmt))
  cat(sprintf("elapsed_s: %.3f\n", t1[["elapsed"]]))

  t2 <- system.time({
    r2 <- extract_sequences(idx, ids_sorted_dense)
  })
  cat(sprintf("case: %s/sorted_unique_dense\n", fmt))
  cat(sprintf("elapsed_s: %.3f\n", t2[["elapsed"]]))

  t3 <- system.time({
    r3 <- extract_sequences(idx, ids_sorted_sparse)
  })
  cat(sprintf("case: %s/sorted_unique_sparse\n", fmt))
  cat(sprintf("elapsed_s: %.3f\n", t3[["elapsed"]]))

  cat(sprintf("rows_random: %d\n", nrow(r1)))
  cat(sprintf("rows_dense: %d\n", nrow(r2)))
  cat(sprintf("rows_sparse: %d\n", nrow(r3)))
}

run_one("fasta")
run_one("fastq")
