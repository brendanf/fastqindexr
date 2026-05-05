test_that("create_index builds FASTA index for single file", {
  path <- tempfile(fileext = ".fa.gz")
  make_fasta_gz(path)

  idx <- create_index(path, type = "fasta")
  expect_s3_class(idx, "fastqindexr_index")
  expect_identical(idx$format, "fasta")
  expect_equal(idx$n_records, 4)
  expect_length(idx$files, 1)
  expect_equal(unname(idx$file_record_counts), 4)
})

test_that("create_index supports concatenated multi-file indexing", {
  p1 <- tempfile(fileext = ".fa.gz")
  p2 <- tempfile(fileext = ".fa.gz")
  make_fasta_gz(p1)
  write_gz_lines(p2, c(">seq5", "ACAC", ">seq6", "TGTG"))

  idx <- create_index(c(p1, p2), type = "fasta")
  expect_equal(idx$n_records, 6)
  expect_equal(unname(idx$file_record_counts), c(4, 2))
  expect_equal(unname(idx$file_record_offsets), c(0, 4, 6))
})

test_that("create_index auto-detects type", {
  path <- tempfile(fileext = ".fq.gz")
  make_fastq_gz(path)
  idx <- create_index(path, type = "auto")
  expect_identical(idx$format, "fastq")
})

test_that("create_index indexes plain FASTA and extract_sequences matches gzip", {
  plain <- tempfile(fileext = ".fa")
  gz <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(c(plain, gz)), add = TRUE)
  writeLines(c(">seq1", "AAAA", ">seq2", "CCCC", ">seq3", "GGGG", ">seq4", "TTTT"), plain)
  write_gz_lines(gz, c(">seq1", "AAAA", ">seq2", "CCCC", ">seq3", "GGGG", ">seq4", "TTTT"))

  idx_plain <- create_index(plain, type = "fasta")
  idx_gz <- create_index(gz, type = "fasta")
  expect_identical(idx_plain$file_compression, "plain")
  expect_identical(idx_gz$file_compression, "gzip")

  ids <- c(3L, 1L, 4L)
  out_plain <- extract_sequences(idx_plain, ids, mode = "indexed")
  out_gz <- extract_sequences(idx_gz, ids, mode = "indexed")
  expect_equal(out_plain, out_gz)
})
