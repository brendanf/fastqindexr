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
