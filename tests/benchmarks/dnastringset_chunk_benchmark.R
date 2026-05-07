#!/usr/bin/env Rscript

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Install Biostrings to run this benchmark.")
}

pkgload::load_all(".", quiet = TRUE)

source(file.path("tests", "benchmarks", "benchmark_sink.R"))
bench_sink_start("dnastringset_chunk_benchmark")
on.exit(bench_sink_stop(), add = TRUE)


estimate_mem_mb <- function(expr) {
  gc(reset = TRUE)
  before <- sum(gc()[, 2L])
  value <- force(expr)
  after <- sum(gc()[, 2L])
  list(value = value, mem_mb = max(0, after - before))
}

extract_via_tempfile <- function(index, ids) {
  tmp <- tempfile(fileext = ".fa")
  on.exit(unlink(tmp), add = TRUE)
  extract_sequences_to_file(
    index,
    ids,
    outfile = tmp,
    type = "fasta",
    append = FALSE
  )
  Biostrings::readDNAStringSet(tmp)
}

run_case <- function(case_name, idx, ids) {
  cat(sprintf("\n[benchmark]\n%s\n", case_name))
  cat(sprintf("records: %d\n", length(ids)))

  t_chunk <- system.time({
    chunk <- estimate_mem_mb(extract_sequences_dnastringset(
      idx,
      ids,
      renumber = "none"
    ))
  })
  cat("method: chunked\n")
  cat(sprintf("elapsed_s: %.3f\n", t_chunk[["elapsed"]]))
  cat(sprintf("mem_delta_mb: %.1f\n", chunk$mem_mb))
  cat(sprintf("n_records_out: %d\n", length(chunk$value)))

  t_tmp <- system.time({
    tmp <- estimate_mem_mb(extract_via_tempfile(idx, ids))
  })
  cat("method: tempfile\n")
  cat(sprintf("elapsed_s: %.3f\n", t_tmp[["elapsed"]]))
  cat(sprintf("mem_delta_mb: %.1f\n", tmp$mem_mb))
  cat(sprintf("n_records_out: %d\n", length(tmp$value)))

  stopifnot(identical(as.character(chunk$value), as.character(tmp$value)))
  invisible(NULL)
}

set.seed(1L)
tmp1 <- tempfile(fileext = ".fa.gz")
tmp2 <- tempfile(fileext = ".fa.gz")
on.exit(unlink(c(tmp1, tmp2)), add = TRUE)
make_benchmark_fasta(tmp1, n = 1200000L, width = 80L)
make_benchmark_fasta(tmp2, n = 80000L, width = 120L)
idx_single <- create_index(tmp1, type = "fasta")
idx_multi <- create_index(c(tmp1, tmp2), type = "fasta")

ids_ordered <- seq_len(1000000L)
ids_random_dup <- sample.int(1200000L, 1000000L, replace = TRUE)
ids_cross_boundary <- c(119990:120000, 120001:120050, 200000, 199999, 120001)

run_case("FASTA/single file ordered", idx_single, ids_ordered)
run_case("FASTA/single file random with duplicates", idx_single, ids_random_dup)
run_case("FASTA/multi-file boundary crossing", idx_multi, ids_cross_boundary)

cat("\n[benchmark]\nFASTA/chunk-size sensitivity\n")
idx_sense <- create_index(tmp1, type = "fasta")
for (target in c(1e6, 5e6, 1e7, 5e7, 1e8)) {
  out <- tryCatch(
    {
      cat(sprintf("chunk size=%g chars n=%d\n", target, length(ids_ordered)))
      t <- system.time({
        x <- extract_sequences_dnastringset(
          idx_sense,
          ids_ordered,
          renumber = "none",
          chunk_chars = target
        )
      })
      if (length(x) != length(ids_ordered)) {
        stop(
          "Extraction returned ",
          length(x),
          " instead of ",
          length(ids_ordered),
          " records.",
          call. = FALSE
        )
      }
      cat(sprintf("elapsed=%.3fs\n", t[["elapsed"]]))
    },
    error = function(e) {
      cat(sprintf("target=%g chars ERROR: %s\n", target, conditionMessage(e)))
    }
  )
}

tmpq1 <- tempfile(fileext = ".fq.gz")
tmpq2 <- tempfile(fileext = ".fq.gz")
on.exit(unlink(c(tmpq1, tmpq2)), add = TRUE)
make_benchmark_fastq(tmpq1, n = 1200000L, width = 80L)
make_benchmark_fastq(tmpq2, n = 80000L, width = 120L)
idxq_single <- create_index(tmpq1, type = "fastq")
idxq_multi <- create_index(c(tmpq1, tmpq2), type = "fastq")

run_case("FASTQ/single file ordered", idxq_single, ids_ordered)
run_case(
  "FASTQ/single file random with duplicates",
  idxq_single,
  ids_random_dup
)
run_case("FASTQ/multi-file boundary crossing", idxq_multi, ids_cross_boundary)
