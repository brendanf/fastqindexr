suppressPackageStartupMessages({
  library(fastqindexr)
})

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Install Biostrings to run this benchmark.")
}

set.seed(1)

tmp <- tempfile(fileext = ".fa.gz")
make_benchmark_fasta(tmp, n = 10000, width = 120)

ids <- sample.int(10000, 2000, replace = TRUE)

idx <- create_index(tmp, type = "fasta")

cat("fastqindexr out-of-order sparse extraction\n")
print(system.time({
  res_fastqindexr <- extract_sequences(idx, ids)
}))

cat("\nBiostrings fasta.index baseline + indexed extraction\n")
print(system.time({
  bi_index <- Biostrings::fasta.index(tmp, seqtype = "DNA")
  selected <- bi_index[ids, , drop = FALSE]
  res_biostrings <- Biostrings::readDNAStringSet(selected)
}))

cat("\nResult cardinalities\n")
cat("fastqindexr:", nrow(res_fastqindexr), "\n")
cat("Biostrings :", length(res_biostrings), "\n")
