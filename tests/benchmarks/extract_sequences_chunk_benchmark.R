#!/usr/bin/env Rscript
# Benchmark `extract_sequences()` after the C++ backend processes unique record
# IDs in fixed-size chunks (see `kChunkUniqueIds` in `cpp_extract_sequences`).
# This is separate from `extract_sequences_dnastringset(..., chunk_chars = )`.
# That argument is an R-visible tuning knob: it must be a finite positive scalar
# (`NA` and `Inf` are errors). Run from package root:
# Rscript tests/benchmarks/extract_sequences_chunk_benchmark.R
#
# After the main timings, runs a **small matrix** of
# `fastqindexr.max_bridge_gap` x `fastqindexr.max_region_bytes` on **moderate**
# request sizes (not full 300k random/permuted workloads), so values stay in a
# reasonable range without pathological tiny caps.

pkgload::load_all(".", quiet = TRUE)

bench_case <- function(label, idx, ids, return = "data.frame") {
  gc(reset = TRUE)
  t <- system.time({
    out <- extract_sequences(idx, ids, return = return)
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
  cat(sprintf(
    "%s return=%s elapsed=%.3fs n=%d%s\n",
    label,
    return,
    t[["elapsed"]],
    length(ids),
    diag_txt
  ))
  invisible(out)
}

bench_with_tuning <- function(
  label,
  idx,
  ids,
  return = "data.frame",
  max_bridge_gap = 64,
  max_region_bytes = 2147483647,
  extract_mode = "indexed",
  diagnostics = TRUE
) {
  old_opts <- options(
    fastqindexr.max_bridge_gap = max_bridge_gap,
    fastqindexr.max_region_bytes = max_region_bytes,
    fastqindexr.extract_mode = extract_mode,
    fastqindexr.extract_diagnostics = diagnostics
  )
  on.exit(options(old_opts), add = TRUE)
  bench_case(label, idx, ids, return = return)
}

bench_region_matrix <- function(idx, label_prefix, ids, tuning_grid) {
  for (i in seq_len(nrow(tuning_grid))) {
    bridge <- tuning_grid$max_bridge_gap[[i]]
    region <- tuning_grid$max_region_bytes[[i]]
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
    extract_mode = "indexed"
  )
  bench_with_tuning(
    sprintf("%s mode=sequential_only", label_prefix),
    idx,
    ids,
    return = "data.frame",
    extract_mode = "sequential_only"
  )
}

cat("generating test sequences\n")
set.seed(1L)
tmp <- tempfile(fileext = ".fa.gz")
on.exit(unlink(tmp), add = TRUE)
make_benchmark_fasta(tmp, n = 1000000L, width = 600L)
cat("creating index\n")
idx <- create_index(tmp, type = "fasta")

cat("extract_sequences chunked-unique backend benchmark\n")

ordered_unique <- seq_len(300000L)
bench_with_tuning("ordered_unique/default", idx, ordered_unique, "data.frame")
bench_with_tuning("ordered_unique/default", idx, ordered_unique, "seq")

cat("execution-mode baseline comparison\n")
set.seed(21L)
bench_mode_compare(
  idx,
  "mode_compare_perm10",
  sample.int(500000L, 10L, replace = FALSE)
)
set.seed(22L)
bench_mode_compare(
  idx,
  "mode_compare_perm1000",
  sample.int(500000L, 1000L, replace = FALSE)
)

tuning_grid <- expand.grid(
  max_bridge_gap = c(65536, 131072, 262144, 524288, 1048576),
  max_region_bytes = c(2000000, 2147483647),
  stringsAsFactors = FALSE
)

cat("region-merge matrix (bridge x max_region_bytes)\n")

set.seed(13L)
bench_region_matrix(
  idx,
  "matrix_perm100k",
  sample.int(500000L, 100000L, replace = FALSE),
  tuning_grid
)
set.seed(13L)
bench_region_matrix(
  idx,
  "matrix_perm10k",
  sample.int(500000L, 10000L, replace = FALSE),
  tuning_grid
)
set.seed(13L)
bench_region_matrix(
  idx,
  "matrix_perm1k",
  sample.int(500000L, 1000L, replace = FALSE),
  tuning_grid
)

set.seed(13L)
bench_region_matrix(
  idx,
  "matrix_perm100",
  sample.int(500000L, 100L, replace = FALSE),
  tuning_grid
)

set.seed(13L)
bench_region_matrix(
  idx,
  "matrix_perm10",
  sample.int(500000L, 10L, replace = FALSE),
  tuning_grid
)

set.seed(2L)
random_dup <- sample.int(500000L, 300000L, replace = TRUE)
bench_with_tuning(
  "random_with_replacement/default",
  idx,
  random_dup,
  "data.frame"
)
bench_with_tuning("random_with_replacement/default", idx, random_dup, "seq")

set.seed(3L)
perm_unique <- sample.int(500000L, 300000L, replace = FALSE)
bench_with_tuning("permuted_unique/default", idx, perm_unique, "data.frame")
bench_with_tuning("permuted_unique/default", idx, perm_unique, "seq")
