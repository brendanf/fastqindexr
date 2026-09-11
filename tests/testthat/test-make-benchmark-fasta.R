test_that("make_benchmark_fasta writes plain FASTA when path is not .gz", {
  path <- tempfile(fileext = ".fa")
  on.exit(unlink(path), add = TRUE)

  make_benchmark_fasta(path, n = 3L, width = 5L)

  expect_true(file.exists(path))
  lines <- readLines(path, warn = FALSE)
  expect_length(lines, 6L)
  expect_true(all(startsWith(lines[c(1L, 3L, 5L)], ">seq")))
})

test_that("make_benchmark_fasta writes gzip FASTA when path ends in .gz", {
  path <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(path), add = TRUE)

  make_benchmark_fasta(path, n = 3L, width = 5L)

  expect_true(file.exists(path))
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  expect_length(lines, 6L)
  expect_true(all(startsWith(lines[c(1L, 3L, 5L)], ">seq")))
})

test_that("make_benchmark_fasta validates path, n, width, alphabet", {
  expect_error(make_benchmark_fasta(c("a", "b")), "`path` must be")
  expect_error(make_benchmark_fasta(NA_character_), "`path` must be")
  expect_error(make_benchmark_fasta(""), "`path` must be")

  ok_path <- tempfile(fileext = ".fa")
  on.exit(unlink(ok_path), add = TRUE)

  expect_error(make_benchmark_fasta(ok_path, n = 0), "`n` must be")
  expect_error(make_benchmark_fasta(ok_path, n = NA_real_), "`n` must be")
  expect_error(
    make_benchmark_fasta(ok_path, n = 1, width = 0),
    "`width` must be"
  )
  expect_error(
    make_benchmark_fasta(ok_path, n = 1, width = 1, alphabet = character()),
    "`alphabet` must contain"
  )
  expect_error(
    make_benchmark_fasta(ok_path, n = 1, width = 1, alphabet = c("A", NA)),
    "`alphabet` must contain"
  )
  expect_error(
    make_benchmark_fasta(ok_path, n = 1, width = 1, alphabet = c("A", "")),
    "`alphabet` must contain"
  )
})

test_that("make_benchmark_fastq writes plain FASTQ when path is not .gz", {
  path <- tempfile(fileext = ".fq")
  on.exit(unlink(path), add = TRUE)

  make_benchmark_fastq(path, n = 3L, width = 5L)

  expect_true(file.exists(path))
  lines <- readLines(path, warn = FALSE)
  expect_length(lines, 12L)
  expect_true(all(startsWith(lines[seq(1L, 12L, by = 4L)], "@seq")))
  expect_true(all(lines[seq(3L, 12L, by = 4L)] == "+"))
  expect_true(all(nchar(lines[seq(2L, 12L, by = 4L)]) == 5L))
  expect_true(all(nchar(lines[seq(4L, 12L, by = 4L)]) == 5L))
})

test_that("make_benchmark_fastq writes gzip FASTQ when path ends in .gz", {
  path <- tempfile(fileext = ".fq.gz")
  on.exit(unlink(path), add = TRUE)

  make_benchmark_fastq(path, n = 2L, width = 4L)

  expect_true(file.exists(path))
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  expect_length(lines, 8L)
  expect_equal(lines[c(1L, 5L)], c("@seq00001", "@seq00002"))
  expect_equal(lines[c(3L, 7L)], c("+", "+"))
})

test_that("make_benchmark_fastq output is round-trippable via create_index", {
  path <- tempfile(fileext = ".fq.gz")
  on.exit(unlink(path), add = TRUE)
  make_benchmark_fastq(path, n = 4L, width = 6L)

  idx <- create_index(path, type = "fastq")
  expect_equal(idx$n_records, 4)
  out <- extract_sequences(idx, c(3L, 1L))
  expect_equal(out$seq_id, c("seq00003", "seq00001"))
  expect_true(all(nchar(out$seq) == 6L))
  expect_true(all(nchar(out$qual) == 6L))
})

test_that("make_benchmark_fastq validates path, n, width, alphabet, qual_alphabet", {
  expect_error(make_benchmark_fastq(c("a", "b")), "`path` must be")
  expect_error(make_benchmark_fastq(NA_character_), "`path` must be")
  expect_error(make_benchmark_fastq(""), "`path` must be")

  ok_path <- tempfile(fileext = ".fq")
  on.exit(unlink(ok_path), add = TRUE)

  expect_error(make_benchmark_fastq(ok_path, n = 0), "`n` must be")
  expect_error(make_benchmark_fastq(ok_path, n = NA_real_), "`n` must be")
  expect_error(
    make_benchmark_fastq(ok_path, n = 1, width = 0),
    "`width` must be"
  )
  expect_error(
    make_benchmark_fastq(ok_path, n = 1, width = 1, alphabet = character()),
    "`alphabet` must contain"
  )
  expect_error(
    make_benchmark_fastq(ok_path, n = 1, width = 1, alphabet = c("A", NA)),
    "`alphabet` must contain"
  )
  expect_error(
    make_benchmark_fastq(ok_path, n = 1, width = 1, alphabet = c("A", "")),
    "`alphabet` must contain"
  )
  expect_error(
    make_benchmark_fastq(
      ok_path,
      n = 1,
      width = 1,
      qual_alphabet = character()
    ),
    "`qual_alphabet` must contain"
  )
  expect_error(
    make_benchmark_fastq(
      ok_path,
      n = 1,
      width = 1,
      qual_alphabet = c("!", NA)
    ),
    "`qual_alphabet` must contain"
  )
  expect_error(
    make_benchmark_fastq(
      ok_path,
      n = 1,
      width = 1,
      qual_alphabet = c("!", "")
    ),
    "`qual_alphabet` must contain"
  )
})
