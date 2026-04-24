test_that("extract_sequences_to_file writes FASTA and preserves order", {
  in_path <- tempfile(fileext = ".fa.gz")
  out_path <- tempfile(fileext = ".fa")
  make_fasta_gz(in_path)
  idx <- create_index(in_path, type = "fasta")

  returned <- extract_sequences_to_file(idx, c(3, 1, 3), outfile = out_path)

  expect_identical(returned, normalizePath(out_path, winslash = "/", mustWork = FALSE))
  expect_equal(
    readLines(out_path, warn = FALSE),
    c(">seq3", "GGGG", ">seq1", "AAAA", ">seq3", "GGGG")
  )
})

test_that("extract_sequences_to_file writes FASTQ and supports gzip output", {
  in_path <- tempfile(fileext = ".fq.gz")
  out_path <- tempfile(fileext = ".fq.gz")
  make_fastq_gz(in_path)
  idx <- create_index(in_path, type = "fastq")

  extract_sequences_to_file(idx, c(4, 2), outfile = out_path, type = "auto")
  con <- gzfile(out_path, "rt")
  on.exit(close(con), add = TRUE)

  expect_equal(
    readLines(con, warn = FALSE),
    c("@r4", "NANA", "+", "%%%%", "@r2", "TTAA", "+", "####")
  )
})

test_that("extract_sequences_to_file supports append and overwrite", {
  in_path <- tempfile(fileext = ".fa.gz")
  out_path <- tempfile(fileext = ".fa")
  make_fasta_gz(in_path)
  idx <- create_index(in_path, type = "fasta")

  extract_sequences_to_file(idx, c(1), outfile = out_path, append = FALSE)
  extract_sequences_to_file(idx, c(2), outfile = out_path, append = TRUE)
  expect_equal(readLines(out_path, warn = FALSE), c(">seq1", "AAAA", ">seq2", "CCCC"))

  extract_sequences_to_file(idx, c(4), outfile = out_path, append = FALSE)
  expect_equal(readLines(out_path, warn = FALSE), c(">seq4", "TTTT"))
})

test_that("extract_sequences_to_file supports .fqi and multifile paths", {
  fqi <- testthat::test_path("fixtures", "cli_fixture.fastq.gz.fqi")
  fq <- testthat::test_path("fixtures", "cli_fixture.fastq.gz")
  out_one <- tempfile(fileext = ".fa")

  extract_sequences_to_file(fqi, c(4, 1), file = fq, outfile = out_one, type = "fasta")
  expect_equal(readLines(out_one, warn = FALSE), c(">r4", "CCCC", ">r1", "ACGT"))

  part1 <- tempfile(fileext = ".part1.fastq.gz")
  part2 <- tempfile(fileext = ".part2.fastq.gz")
  make_fastq_gz(part1)
  write_gz_lines(part2, c("@x1", "CCCC", "+", "!!!!", "@x2", "GGGG", "+", "####"))
  idx <- create_index(c(part1, part2), type = "fastq")
  out_two <- tempfile(fileext = ".fq")

  extract_sequences_to_file(idx, c(5, 2, 6), outfile = out_two, type = "fastq")
  expect_equal(
    readLines(out_two, warn = FALSE),
    c("@x1", "CCCC", "+", "!!!!", "@r2", "TTAA", "+", "####", "@x2", "GGGG", "+", "####")
  )
})

test_that("extract_sequences_to_file strict-increasing output matches in-memory", {
  tmp_in <- tempfile(fileext = ".fa.gz")
  out_stream <- tempfile(fileext = ".fa")
  out_ref <- tempfile(fileext = ".fa")
  on.exit(unlink(c(tmp_in, out_stream, out_ref)), add = TRUE)
  make_benchmark_fasta(tmp_in, n = 200L, width = 30L)
  idx <- create_index(tmp_in, type = "fasta")
  set.seed(7L)
  ids <- sort(sample.int(200L, 80L, replace = FALSE))

  extract_sequences_to_file(idx, ids, outfile = out_stream, append = FALSE)

  ex <- extract_sequences(idx, ids)
  ref_lines <- c(rbind(paste0(">", ex$seq_id), ex$seq))
  writeLines(ref_lines, out_ref)
  expect_equal(readLines(out_stream, warn = FALSE), readLines(out_ref, warn = FALSE))
})
