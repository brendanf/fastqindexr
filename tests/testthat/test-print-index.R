test_that("print.fastqindexr_index summarizes format, files, and records", {
  path <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(path), add = TRUE)
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")

  out <- capture.output(res <- print(idx))
  expect_identical(res, idx)
  expect_match(out, "<fastqindexr_index>", all = FALSE)
  expect_match(out, "format: fasta", all = FALSE)
  expect_match(out, "files: 1", all = FALSE)
  expect_match(out, "records: 4", all = FALSE)
})

test_that("print.fastqindexr_index reports multi-file FASTQ index", {
  p1 <- tempfile(fileext = ".fq.gz")
  p2 <- tempfile(fileext = ".fq.gz")
  on.exit(unlink(c(p1, p2)), add = TRUE)
  make_fastq_gz(p1)
  write_gz_lines(p2, c("@x1", "AAAA", "+", "!!!!"))
  idx <- create_index(c(p1, p2), type = "fastq")

  out <- capture.output(print(idx))
  expect_match(out, "format: fastq", all = FALSE)
  expect_match(out, "files: 2", all = FALSE)
  expect_match(out, "records: 5", all = FALSE)
})
