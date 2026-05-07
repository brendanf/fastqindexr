# Optional benchmark: direct extract-to-file vs extract-to-list + writeLines.
# Not run by testthat.
#
# Usage (from package root, so current sources can be loaded):
#   Rscript tests/benchmarks/extract_to_file_vs_list_benchmark.R
# If the pkgload package is available, the script loads the in-tree build via
# load_all(); otherwise use a recent `install` of `fastqindexr`.
if (file.exists("DESCRIPTION") && file.exists("R") && file.exists("src") &&
    requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(fastqindexr))
}

source(file.path("tests", "benchmarks", "benchmark_sink.R"))
bench_sink_start("extract_to_file_vs_list_benchmark")
on.exit(bench_sink_stop(), add = TRUE)

# Flatten a FASTA list (seq_id, seq) to lines like extract_sequences_to_file
# for plain FASTA output.
write_fasta_list_to_file <- function(L, path) {
  lines <- c(rbind(paste0(">", L$seq_id), L$seq))
  writeLines(lines, path)
}

write_fastq_list_to_file <- function(L, path) {
  n <- length(L$seq)
  lines <- character(n * 4L)
  for (i in seq_len(n)) {
    j <- (i - 1L) * 4L + 1L
    lines[j] <- paste0("@", L$seq_id[[i]])
    lines[j + 1L] <- L$seq[[i]]
    lines[j + 2L] <- "+"
    lines[j + 3L] <- L$qual[[i]]
  }
  writeLines(lines, path)
}

set.seed(1L)
n_rec <- 10000L
k <- 2000L
tmp_in <- tempfile(fileext = ".fa.gz")
out_direct <- tempfile(fileext = ".fa")
out_via_list <- tempfile(fileext = ".fa")
on.exit(
  unlink(c(tmp_in, out_direct, out_via_list)),
  add = TRUE
)

run_one <- function(fmt) {
  if (fmt == "fasta") {
    make_benchmark_fasta(tmp_in, n = n_rec, width = 120L)
    out_direct <<- tempfile(fileext = ".fa")
    out_via_list <<- tempfile(fileext = ".fa")
    out_type <- "fasta"
  } else {
    make_benchmark_fastq(tmp_in, n = n_rec, width = 120L)
    out_direct <<- tempfile(fileext = ".fq")
    out_via_list <<- tempfile(fileext = ".fq")
    out_type <- "fastq"
  }

  ids <- sort(sample.int(n_rec, k, replace = FALSE))
  cat(sprintf("\n[benchmark]\n%s: build index\n", toupper(fmt)))
  t_idx <- system.time({
    idx <- create_index(tmp_in, type = fmt)
  })
  cat(sprintf("case: %s/index\n", fmt))
  cat(sprintf("elapsed_s: %.3f\n", t_idx[["elapsed"]]))

  t_direct <- system.time({
    extract_sequences_to_file(
      idx,
      seq_idx = ids,
      outfile = out_direct,
      type = out_type,
      append = FALSE,
      compress = FALSE
    )
  })
  cat(sprintf("case: %s/direct_to_file\n", fmt))
  cat(sprintf("elapsed_s: %.3f\n", t_direct[["elapsed"]]))

  t_list <- system.time({
    L <- extract_sequences(idx, seq_idx = ids, return = "list")
    if (fmt == "fasta") {
      write_fasta_list_to_file(L, out_via_list)
    } else {
      write_fastq_list_to_file(L, out_via_list)
    }
  })
  cat(sprintf("case: %s/list_then_write\n", fmt))
  cat(sprintf("elapsed_s: %.3f\n", t_list[["elapsed"]]))

  cat(sprintf("size_direct_bytes: %d\n", file.size(out_direct)))
  cat(sprintf("size_list_bytes: %d\n", file.size(out_via_list)))
  if (
    !identical(
      readLines(out_direct, warn = FALSE),
      readLines(out_via_list, warn = FALSE)
    )
  ) {
    stop("Output mismatch between paths (benchmark sanity check).")
  }
  cat("outputs_match: TRUE\n")
}

run_one("fasta")
run_one("fastq")
