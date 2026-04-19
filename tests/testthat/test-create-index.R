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
