#!/usr/bin/env Rscript
# After default partition timings, runs a **small matrix** of
# `fastqindexr.max_bridge_gap` x `fastqindexr.max_region_records` on the same
# partitions (moderate grid, no pathological tiny caps).

pkgload::load_all(".", quiet = TRUE)

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
  elapsed <- system.time({
    ret <- extract_sequences_to_file(
      idx,
      partitions,
      outfile = outs,
      append = FALSE
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
  cat(sprintf("%s elapsed: %.3fs%s\n", label, elapsed, diag_txt))
}

bench_partition_with_tuning <- function(
  label,
  partitions,
  max_bridge_gap = 64,
  max_region_records = 2147483647,
  extract_mode = "indexed",
  diagnostics = TRUE
) {
  old_opts <- options(
    fastqindexr.max_bridge_gap = max_bridge_gap,
    fastqindexr.max_region_records = max_region_records,
    fastqindexr.extract_mode = extract_mode,
    fastqindexr.extract_diagnostics = diagnostics
  )
  on.exit(options(old_opts), add = TRUE)
  bench_partition(label, partitions)
}

bench_partition_matrix <- function(label_prefix, partitions, tuning_grid) {
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
    bench_partition_with_tuning(lab, partitions, bridge, region)
  }
}

bench_partition_mode_compare <- function(label_prefix, partitions) {
  bench_partition_with_tuning(
    sprintf("%s mode=indexed", label_prefix),
    partitions,
    extract_mode = "indexed"
  )
  bench_partition_with_tuning(
    sprintf("%s mode=sequential_only", label_prefix),
    partitions,
    extract_mode = "sequential_only"
  )
}

cat("Partitioned strict-increasing extraction benchmark\n")
bench_partition_with_tuning("contiguous-8/default", cont)
bench_partition_with_tuning("round_robin-8/default", rr)
cat("execution-mode baseline comparison\n")
bench_partition_mode_compare("mode_compare_contiguous-8", cont)
bench_partition_mode_compare("mode_compare_round_robin-8", rr)

tuning_grid <- expand.grid(
  max_bridge_gap = c(0, 64, 256),
  max_region_records = c(500000, 2000000, 2147483647),
  stringsAsFactors = FALSE
)

cat("region-merge matrix (partitioned contiguous vs round_robin)\n")
bench_partition_matrix("matrix_contig-8", cont, tuning_grid)
bench_partition_matrix("matrix_rr-8", rr, tuning_grid)
