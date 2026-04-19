write_gz_lines <- function(path, lines) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con)
}

make_fasta_gz <- function(path) {
  write_gz_lines(path, c(
    ">seq1", "AAAA",
    ">seq2", "CCCC",
    ">seq3", "GGGG",
    ">seq4", "TTTT"
  ))
}

make_fastq_gz <- function(path) {
  write_gz_lines(path, c(
    "@r1", "ACGT", "+", "!!!!",
    "@r2", "TTAA", "+", "####",
    "@r3", "GCGC", "+", "$$$$",
    "@r4", "NANA", "+", "%%%%"
  ))
}
