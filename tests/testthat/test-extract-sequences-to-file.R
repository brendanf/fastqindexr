test_that("extract_sequences_to_file writes FASTA and preserves order", {
  in_path <- tempfile(fileext = ".fa.gz")
  out_path <- tempfile(fileext = ".fa")
  make_fasta_gz(in_path)
  idx <- create_index(in_path, type = "fasta")

  returned <- extract_sequences_to_file(idx, c(3, 1, 3), outfile = out_path)

  expect_identical(
    returned,
    normalizePath(out_path, winslash = "/", mustWork = FALSE)
  )
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
  expect_equal(
    readLines(out_path, warn = FALSE),
    c(">seq1", "AAAA", ">seq2", "CCCC")
  )

  extract_sequences_to_file(idx, c(4), outfile = out_path, append = FALSE)
  expect_equal(readLines(out_path, warn = FALSE), c(">seq4", "TTTT"))
})

test_that("extract_sequences_to_file supports .fqi and multifile paths", {
  fqi <- testthat::test_path("fixtures", "cli_fixture.fastq.gz.fqi")
  fq <- testthat::test_path("fixtures", "cli_fixture.fastq.gz")
  out_one <- tempfile(fileext = ".fa")

  extract_sequences_to_file(
    fqi,
    c(4, 1),
    file = fq,
    outfile = out_one,
    type = "fasta"
  )
  expect_equal(
    readLines(out_one, warn = FALSE),
    c(">r4", "CCCC", ">r1", "ACGT")
  )

  part1 <- tempfile(fileext = ".part1.fastq.gz")
  part2 <- tempfile(fileext = ".part2.fastq.gz")
  make_fastq_gz(part1)
  write_gz_lines(
    part2,
    c("@x1", "CCCC", "+", "!!!!", "@x2", "GGGG", "+", "####")
  )
  idx <- create_index(c(part1, part2), type = "fastq")
  out_two <- tempfile(fileext = ".fq")

  extract_sequences_to_file(idx, c(5, 2, 6), outfile = out_two, type = "fastq")
  expect_equal(
    readLines(out_two, warn = FALSE),
    c(
      "@x1",
      "CCCC",
      "+",
      "!!!!",
      "@r2",
      "TTAA",
      "+",
      "####",
      "@x2",
      "GGGG",
      "+",
      "####"
    )
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
  expect_equal(
    readLines(out_stream, warn = FALSE),
    readLines(out_ref, warn = FALSE)
  )
})

test_that("extract_sequences_to_file supports partitioned seq_idx and outfiles", {
  in_path <- tempfile(fileext = ".fa.gz")
  out1 <- tempfile(fileext = ".fa")
  out2 <- tempfile(fileext = ".fa")
  make_fasta_gz(in_path)
  idx <- create_index(in_path, type = "fasta")

  returned <- extract_sequences_to_file(
    idx,
    seq_idx = list(c(3, 1, 3), c(2, 4)),
    outfile = c(out1, out2)
  )
  expect_equal(
    returned,
    normalizePath(c(out1, out2), winslash = "/", mustWork = FALSE)
  )
  expect_equal(
    readLines(out1, warn = FALSE),
    c(">seq3", "GGGG", ">seq1", "AAAA", ">seq3", "GGGG")
  )
  expect_equal(
    readLines(out2, warn = FALSE),
    c(">seq2", "CCCC", ">seq4", "TTTT")
  )
})

test_that("extract_sequences_to_file partitioned mode handles empty partitions", {
  in_path <- tempfile(fileext = ".fa.gz")
  out1 <- tempfile(fileext = ".fa")
  out2 <- tempfile(fileext = ".fa")
  make_fasta_gz(in_path)
  idx <- create_index(in_path, type = "fasta")

  writeLines("stale", out1)
  writeLines("keep", out2)
  extract_sequences_to_file(
    idx,
    seq_idx = list(integer(), integer()),
    outfile = c(out1, out2),
    append = FALSE
  )
  expect_equal(readLines(out1, warn = FALSE), character())
  expect_equal(readLines(out2, warn = FALSE), character())

  writeLines("keep1", out1)
  writeLines("keep2", out2)
  extract_sequences_to_file(
    idx,
    seq_idx = list(integer(), integer()),
    outfile = c(out1, out2),
    append = TRUE
  )
  expect_equal(readLines(out1, warn = FALSE), "keep1")
  expect_equal(readLines(out2, warn = FALSE), "keep2")
})

test_that("extract_sequences_to_file partitioned mode supports compress vector", {
  in_path <- tempfile(fileext = ".fq.gz")
  out1 <- tempfile(fileext = ".fq.gz")
  out2 <- tempfile(fileext = ".fq")
  make_fastq_gz(in_path)
  idx <- create_index(in_path, type = "fastq")

  extract_sequences_to_file(
    idx,
    seq_idx = list(c(1, 2), c(4)),
    outfile = c(out1, out2),
    type = "fastq",
    compress = c(TRUE, FALSE)
  )
  con <- gzfile(out1, "rt")
  on.exit(close(con), add = TRUE)
  expect_equal(
    readLines(con, warn = FALSE),
    c("@r1", "ACGT", "+", "!!!!", "@r2", "TTAA", "+", "####")
  )
  expect_equal(readLines(out2, warn = FALSE), c("@r4", "NANA", "+", "%%%%"))
})

test_that("partitioned extraction matches in-memory for contiguous and round_robin", {
  tmp_in <- tempfile(fileext = ".fa.gz")
  out_c1 <- tempfile(fileext = ".fa")
  out_c2 <- tempfile(fileext = ".fa")
  out_r1 <- tempfile(fileext = ".fa")
  out_r2 <- tempfile(fileext = ".fa")
  on.exit(unlink(c(tmp_in, out_c1, out_c2, out_r1, out_r2)), add = TRUE)
  make_benchmark_fasta(tmp_in, n = 200L, width = 30L)
  idx <- create_index(tmp_in, type = "fasta")
  ids <- sort(sample.int(200L, 80L, replace = FALSE))

  cont <- partition_seq_idx(ids, n_parts = 2, strategy = "contiguous")
  rr <- partition_seq_idx(ids, n_parts = 2, strategy = "round_robin")

  extract_sequences_to_file(
    idx,
    cont,
    outfile = c(out_c1, out_c2),
    append = FALSE
  )
  extract_sequences_to_file(
    idx,
    rr,
    outfile = c(out_r1, out_r2),
    append = FALSE
  )

  ex_c1 <- extract_sequences(idx, cont[[1]])
  ex_c2 <- extract_sequences(idx, cont[[2]])
  ex_r1 <- extract_sequences(idx, rr[[1]])
  ex_r2 <- extract_sequences(idx, rr[[2]])
  expect_equal(
    readLines(out_c1, warn = FALSE),
    c(rbind(paste0(">", ex_c1$seq_id), ex_c1$seq))
  )
  expect_equal(
    readLines(out_c2, warn = FALSE),
    c(rbind(paste0(">", ex_c2$seq_id), ex_c2$seq))
  )
  expect_equal(
    readLines(out_r1, warn = FALSE),
    c(rbind(paste0(">", ex_r1$seq_id), ex_r1$seq))
  )
  expect_equal(
    readLines(out_r2, warn = FALSE),
    c(rbind(paste0(">", ex_r2$seq_id), ex_r2$seq))
  )
})
