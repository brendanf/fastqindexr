#!/usr/bin/env Rscript
# After default partition timings, runs a **small matrix** of
# `fastqindexr.max_bridge_gap` x `fastqindexr.max_region_records` on the same
# partitions (moderate grid, no pathological tiny caps).

pkgload::load_all(".", quiet = TRUE)

source(file.path("tests", "benchmarks", "benchmark_sink.R"))
bench_sink_start("partitioned_extract_to_file_benchmark")
on.exit(bench_sink_stop(), add = TRUE)

bench_cache_file <- function(name) {
  cache_dir <- file.path("tests", "benchmarks", "cache")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  file.path(cache_dir, name)
}

tmp <- bench_cache_file("benchmark_250k_w60.fa.gz")
if (!file.exists(tmp)) {
  cat("generating test sequences (cache miss)\n")
  make_benchmark_fasta(tmp, n = 250000L, width = 60L)
} else {
  cat("using cached test sequences\n")
}
idx <- create_index(tmp, type = "fasta")
idx_gz <- idx

tmp_plain <- bench_cache_file("benchmark_250k_w60.fa")
if (!file.exists(tmp_plain)) {
  cat("creating plain FASTA benchmark cache\n")
  make_benchmark_fasta(tmp_plain, n = 250000L, width = 60L)
}
idx_plain <- create_index(tmp_plain, type = "fasta")

ids <- seq_len(200000L)
cont <- partition_seq_idx(ids, n_parts = 8L, strategy = "contiguous")
rr <- partition_seq_idx(ids, n_parts = 8L, strategy = "round_robin")

bench_partition <- function(label, partitions, mode = "auto", out_ext = ".fa") {
  outs <- replicate(
    length(partitions),
    tempfile(fileext = out_ext),
    simplify = TRUE
  )
  on.exit(unlink(outs), add = TRUE)
  elapsed <- system.time({
    ret <- extract_sequences_to_file(
      idx,
      partitions,
      outfile = outs,
      append = FALSE,
      mode = mode
    )
  })[["elapsed"]]
  diag <- attr(ret, "fastqindexr_diagnostics", exact = TRUE)
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
  cat(sprintf("elapsed_s: %.3f\n", elapsed))
  if (nzchar(diag_txt)) {
    cat(sprintf("diagnostics:%s\n", diag_txt))
  }
}

bench_partition_with_tuning <- function(
  label,
  partitions,
  out_ext = ".fa",
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
  bench_partition(label, partitions, mode = mode, out_ext = out_ext)
}

bench_partition_matrix <- function(
  label_prefix,
  partitions,
  tuning_grid,
  out_ext = ".fa"
) {
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
    bench_partition_with_tuning(lab, partitions, out_ext, bridge, region)
  }
}

bench_partition_mode_compare <- function(
  label_prefix,
  partitions,
  out_ext = ".fa"
) {
  bench_partition_with_tuning(
    sprintf("%s mode=indexed", label_prefix),
    partitions,
    out_ext = out_ext,
    extract_mode = "indexed",
    mode = "auto"
  )
  bench_partition_with_tuning(
    sprintf("%s mode=sequential", label_prefix),
    partitions,
    out_ext = out_ext,
    extract_mode = "indexed",
    mode = "sequential"
  )
}

bench_log_header("FASTA: partitioned strict-increasing extraction benchmark")
bench_partition_with_tuning("fasta/contiguous-8/default", cont)
bench_partition_with_tuning("fasta/round_robin-8/default", rr)
bench_log_header("FASTA: execution-mode baseline comparison")
bench_partition_mode_compare("fasta/mode_compare_contiguous-8", cont)
bench_partition_mode_compare("fasta/mode_compare_round_robin-8", rr)

bench_log_header("FASTA plain: execution-mode baseline comparison")
idx <- idx_plain
bench_partition_mode_compare("fasta_plain/mode_compare_contiguous-8", cont)
bench_partition_mode_compare("fasta_plain/mode_compare_round_robin-8", rr)
idx <- idx_gz

tuning_grid <- expand.grid(
  max_bridge_gap = c(0, 64, 256),
  max_region_records = c(500000, 2000000, 2147483647),
  stringsAsFactors = FALSE
)

bench_log_header("FASTA: region-merge matrix (contiguous vs round_robin)")
bench_partition_matrix("fasta/matrix_contig-8", cont, tuning_grid)
bench_partition_matrix("fasta/matrix_rr-8", rr, tuning_grid)

tmp_q <- bench_cache_file("benchmark_250k_w60.fq.gz")
if (!file.exists(tmp_q)) {
  make_benchmark_fastq(tmp_q, n = 250000L, width = 60L)
}
idx_q <- create_index(tmp_q, type = "fastq")
idx <- idx_q
bench_log_header("FASTQ: partitioned strict-increasing extraction benchmark")
bench_partition_with_tuning("fastq/contiguous-8/default", cont, out_ext = ".fq")
bench_partition_with_tuning("fastq/round_robin-8/default", rr, out_ext = ".fq")
bench_log_header("FASTQ: execution-mode baseline comparison")
bench_partition_mode_compare(
  "fastq/mode_compare_contiguous-8",
  cont,
  out_ext = ".fq"
)
bench_partition_mode_compare(
  "fastq/mode_compare_round_robin-8",
  rr,
  out_ext = ".fq"
)
bench_log_header("FASTQ: region-merge matrix (contiguous vs round_robin)")
bench_partition_matrix(
  "fastq/matrix_contig-8",
  cont,
  tuning_grid,
  out_ext = ".fq"
)
bench_partition_matrix("fastq/matrix_rr-8", rr, tuning_grid, out_ext = ".fq")
