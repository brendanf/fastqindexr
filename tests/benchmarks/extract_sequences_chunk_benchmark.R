#!/usr/bin/env Rscript
# Benchmark `extract_sequences()` with region-planned indexed extraction.
# This is separate from `extract_sequences_dnastringset(..., chunk_chars = )`.
# That argument is an R-visible tuning knob: it must be a finite positive scalar
# (`NA` and `Inf` are errors). Run from package root:
# Rscript tests/benchmarks/extract_sequences_chunk_benchmark.R
#
# After the main timings, runs a **small matrix** of
# `fastqindexr.max_bridge_gap` x `fastqindexr.max_region_records` on **moderate**
# request sizes (not full 300k random/permuted workloads), so values stay in a
# reasonable range without pathological tiny caps.

pkgload::load_all(".", quiet = TRUE)

source(file.path("tests", "benchmarks", "benchmark_sink.R"))
bench_sink_start("extract_sequences_chunk_benchmark")
on.exit(bench_sink_stop(), add = TRUE)

bench_cache_file <- function(name) {
  cache_dir <- file.path("tests", "benchmarks", "cache")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  file.path(cache_dir, name)
}

bench_case <- function(label, idx, ids, return = "data.frame", mode = "auto") {
  gc(reset = TRUE)
  t <- system.time({
    out <- extract_sequences(idx, ids, return = return, mode = mode)
  })
  diag <- attr(out, "fastqindexr_diagnostics", exact = TRUE)
  diag_txt <- ""
  if (!is.null(diag)) {
    diag_txt <- sprintf(
      " regions=%s attempts=%s failures=%s fallback_calls=%s fallback_records=%s",
      diag$regions_planned,
      diag$extract_attempts,
      diag$extract_failures,
      diag$fallback_invocations,
      diag$fallback_records
    )
  }
  cat(sprintf("case: %s\n", label))
  cat(sprintf("    return: %s\n", return))
  cat(sprintf("    elapsed_s: %.3f\n", t[["elapsed"]]))
  cat(sprintf("    n_ids: %d\n", length(ids)))
  if (nzchar(diag_txt)) {
    cat(sprintf("    diagnostics:%s\n", diag_txt))
  }
  invisible(out)
}

bench_with_tuning <- function(
  label,
  idx,
  ids,
  return = "data.frame",
  max_bridge_gap = 64,
  max_region_records = 2147483647,
  extract_mode = "indexed",
  diagnostics = TRUE,
  mode = "auto"
) {
  old_opts <- options(
    fastqindexr.max_bridge_gap = max_bridge_gap,
    fastqindexr.max_region_records = max_region_records,
    fastqindexr.extract_mode = extract_mode,
    fastqindexr.extract_diagnostics = diagnostics
  )
  on.exit(options(old_opts), add = TRUE)
  bench_case(label, idx, ids, return = return, mode = mode)
}

bench_region_matrix <- function(idx, label_prefix, ids, tuning_grid) {
  lab <- sprintf("%s mode=sequential", label_prefix)
  bench_with_tuning(
    lab,
    idx,
    ids,
    "data.frame",
    extract_mode = "indexed",
    mode = "sequential"
  )
  for (i in seq_len(nrow(tuning_grid))) {
    bridge <- tuning_grid$max_bridge_gap[[i]]
    region <- tuning_grid$max_region_records[[i]]
    lab <- sprintf(
      "%s mode=%s b=%s r=%s",
      label_prefix,
      "indexed",
      format(bridge, scientific = FALSE, trim = TRUE),
      format(region, scientific = FALSE, trim = TRUE)
    )
    bench_with_tuning(lab, idx, ids, "data.frame", bridge, region)
  }
}

bench_mode_compare <- function(idx, label_prefix, ids) {
  bench_with_tuning(
    sprintf("%s mode=indexed", label_prefix),
    idx,
    ids,
    return = "data.frame",
    extract_mode = "indexed",
    mode = "auto"
  )
  bench_with_tuning(
    sprintf("%s mode=sequential", label_prefix),
    idx,
    ids,
    return = "data.frame",
    extract_mode = "indexed",
    mode = "sequential"
  )
}

set.seed(1L)
tmp <- bench_cache_file("benchmark_1m_w600.fa.gz")
if (!file.exists(tmp)) {
  cat("generating test sequences (cache miss)\n")
  make_benchmark_fasta(tmp, n = 1000000L, width = 600L)
} else {
  cat("using cached test sequences\n")
}
cat("creating index\n")
idx <- create_index(tmp, type = "fasta")

tmp_plain <- bench_cache_file("benchmark_1m_w600.fa")
if (!file.exists(tmp_plain)) {
  cat("creating plain FASTA benchmark cache\n")
  make_benchmark_fasta(tmp_plain, n = 1000000L, width = 600L)
}
idx_plain <- create_index(tmp_plain, type = "fasta")

bench_log_header("FASTA: extract_sequences chunked backend benchmark")

ordered_unique <- seq_len(300000L)
bench_with_tuning("ordered_unique/default", idx, ordered_unique, "data.frame")
bench_with_tuning("ordered_unique/default", idx, ordered_unique, "seq")

bench_log_header("FASTA: execution-mode baseline comparison")
set.seed(21L)
bench_mode_compare(
  idx,
  "mode_compare_perm10",
  sample.int(1000000L, 10L, replace = FALSE)
)
set.seed(22L)
bench_mode_compare(
  idx,
  "mode_compare_perm1000",
  sample.int(1000000L, 1000L, replace = FALSE)
)

tuning_grid <- expand.grid(
  max_bridge_gap = c(64, 256, 1024, 4096, 16384, 65536),
  max_region_records = c(200000, 2000000),
  stringsAsFactors = FALSE
)

bench_log_header("FASTA: region-merge matrix (bridge x max_region_records)")

set.seed(13L)
bench_region_matrix(
  idx,
  "matrix_perm1M",
  sample.int(1000000L, 1000000L, replace = FALSE),
  tuning_grid
)

set.seed(13L)
bench_region_matrix(
  idx,
  "matrix_perm100k",
  sample.int(1000000L, 100000L, replace = FALSE),
  tuning_grid
)
set.seed(13L)
bench_region_matrix(
  idx,
  "matrix_perm10k",
  sample.int(1000000L, 10000L, replace = FALSE),
  tuning_grid
)
set.seed(13L)
bench_region_matrix(
  idx,
  "matrix_perm1k",
  sample.int(1000000L, 1000L, replace = FALSE),
  tuning_grid
)

set.seed(13L)
bench_region_matrix(
  idx,
  "matrix_perm100",
  sample.int(1000000L, 100L, replace = FALSE),
  tuning_grid
)

set.seed(13L)
bench_region_matrix(
  idx,
  "matrix_perm10",
  sample.int(1000000L, 10L, replace = FALSE),
  tuning_grid
)

set.seed(2L)
random_dup <- sample.int(1000000L, 300000L, replace = TRUE)
bench_with_tuning(
  "random_with_replacement/default",
  idx,
  random_dup,
  "data.frame"
)
bench_with_tuning("random_with_replacement/default", idx, random_dup, "seq")

set.seed(3L)
perm_unique <- sample.int(1000000L, 300000L, replace = FALSE)
bench_with_tuning("permuted_unique/default", idx, perm_unique, "data.frame")
bench_with_tuning("permuted_unique/default", idx, perm_unique, "seq")

bench_log_header("FASTA plain: mode comparison")
set.seed(23L)
bench_mode_compare(
  idx_plain,
  "plain_mode_compare_perm1000",
  sample.int(1000000L, 1000L, replace = FALSE)
)

tmp_q <- bench_cache_file("benchmark_1m_w600.fq.gz")
if (!file.exists(tmp_q)) {
  make_benchmark_fastq(tmp_q, n = 1000000L, width = 600L)
}
idx_q <- create_index(tmp_q, type = "fastq")

bench_log_header("FASTQ: extract_sequences chunked backend benchmark")
bench_with_tuning(
  "fastq/ordered_unique/default",
  idx_q,
  ordered_unique,
  "data.frame"
)
bench_with_tuning("fastq/ordered_unique/default", idx_q, ordered_unique, "seq")

bench_log_header("FASTQ: execution-mode baseline comparison")
set.seed(31L)
bench_mode_compare(
  idx_q,
  "fastq/mode_compare_perm10",
  sample.int(1000000L, 10L, replace = FALSE)
)
set.seed(32L)
bench_mode_compare(
  idx_q,
  "fastq/mode_compare_perm1000",
  sample.int(1000000L, 1000L, replace = FALSE)
)

bench_log_header("FASTQ: region-merge matrix")
set.seed(41L)
bench_region_matrix(
  idx_q,
  "fastq/matrix_perm10k",
  sample.int(1000000L, 10000L, replace = FALSE),
  tuning_grid
)

set.seed(42L)
random_dup_q <- sample.int(1000000L, 300000L, replace = TRUE)
bench_with_tuning(
  "fastq/random_with_replacement/default",
  idx_q,
  random_dup_q,
  "data.frame"
)
bench_with_tuning(
  "fastq/random_with_replacement/default",
  idx_q,
  random_dup_q,
  "seq"
)
