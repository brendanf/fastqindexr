pkgload::load_all(".", quiet = TRUE)

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Install Biostrings to run this benchmark.")
}

source(file.path("tests", "benchmarks", "benchmark_sink.R"))
bench_sink_start("biostrings_fasta_index_benchmark")
on.exit(bench_sink_stop(), add = TRUE)

set.seed(1)

n_records <- 10000L
n_extract <- 2000L

requests <- list(
  contiguous = (n_records %/% 2L) + seq_len(n_extract),
  sorted_unique = sort(sample.int(n_records, n_extract, replace = FALSE)),
  unsorted_dup = sample.int(n_records, n_extract, replace = TRUE)
)

benchmarks <- data.frame(
  compression = character(),
  implementation = character(),
  operation = character(),
  elapsed = numeric(),
  stringsAsFactors = FALSE
)

for (comp in c("gzip", "plain")) {
  ext <- if (comp == "gzip") ".fa.gz" else ".fa"
  tmp <- tempfile(fileext = ext)
  on.exit(unlink(tmp), add = TRUE)
  make_benchmark_fasta(tmp, n = n_records, width = 120L)

  for (impl in c("fastqindexr", "Biostrings")) {
    label_base <- sprintf("%s/%s", comp, impl)
    t_idx <- system.time({
      if (impl == "fastqindexr") {
        idx <- create_index(tmp, type = "fasta")
      } else {
        idx <- Biostrings::fasta.index(tmp, seqtype = "DNA")
      }
    })
    bench_log_timing(sprintf("%s/index", label_base), t_idx)
    benchmarks <- rbind(
      benchmarks,
      data.frame(
        compression = comp,
        implementation = impl,
        operation = "index",
        elapsed = t_idx[["elapsed"]],
        stringsAsFactors = FALSE
      )
    )

    for (op in names(requests)) {
      ids <- requests[[op]]
      t_ex <- system.time({
        if (impl == "fastqindexr") {
          res <- extract_sequences(idx, ids)
          n_out <- nrow(res)
        } else {
          res <- Biostrings::readDNAStringSet(idx[ids, , drop = FALSE])
          n_out <- length(res)
        }
      })
      bench_log_timing(sprintf("%s/%s", label_base, op), t_ex)
      cat(sprintf("n_records_out: %d\n", n_out))
      benchmarks <- rbind(
        benchmarks,
        data.frame(
          compression = comp,
          implementation = impl,
          operation = op,
          elapsed = t_ex[["elapsed"]],
          stringsAsFactors = FALSE
        )
      )
    }
  }
}
