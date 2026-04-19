test_that("extract_sequences preserves input order and duplicates for FASTA", {
  path <- tempfile(fileext = ".fa.gz")
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")

  out <- extract_sequences(idx, c(3, 1, 3))
  expect_identical(names(out), c("seq_id", "seq"))
  expect_equal(out$seq_id, c("seq3", "seq1", "seq3"))
  expect_equal(out$seq, c("GGGG", "AAAA", "GGGG"))
})

test_that("extract_sequences returns quality for FASTQ", {
  path <- tempfile(fileext = ".fq.gz")
  make_fastq_gz(path)
  idx <- create_index(path, type = "fastq")

  out <- extract_sequences(idx, c(4, 2))
  expect_identical(names(out), c("seq_id", "seq", "qual"))
  expect_equal(out$seq_id, c("r4", "r2"))
  expect_equal(out$seq, c("NANA", "TTAA"))
  expect_equal(out$qual, c("%%%%", "####"))
})

test_that("extract_sequences supports file override for moved inputs", {
  old_path <- tempfile(fileext = ".fa.gz")
  new_path <- tempfile(fileext = ".fa.gz")
  make_fasta_gz(old_path)
  idx <- create_index(old_path, type = "fasta")
  file.rename(old_path, new_path)

  out <- extract_sequences(idx, c(2, 4), file = new_path)
  expect_equal(out$seq_id, c("seq2", "seq4"))
})

test_that("extract_sequences validates ids", {
  path <- tempfile(fileext = ".fa.gz")
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")

  expect_error(extract_sequences(idx, c(0, 1)), "positive whole numbers")
  expect_error(extract_sequences(idx, c(1, 99)), "exceed available records")
})
