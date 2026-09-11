test_that("extract_sequences_dnastringset mode='sequential' matches indexed for gzip FASTA", {
  skip_if_not_installed("Biostrings")
  path <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(path), add = TRUE)
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")

  ids <- c(3L, 1L, 3L, 2L)
  ref <- extract_sequences_dnastringset(idx, ids, mode = "indexed")
  got <- extract_sequences_dnastringset(idx, ids, mode = "sequential")
  expect_equal(unname(as.character(got)), unname(as.character(ref)))
  expect_equal(names(got), names(ref))
})

test_that("extract_sequences_dnastringset mode='sequential' works for plain FASTA", {
  skip_if_not_installed("Biostrings")
  path <- tempfile(fileext = ".fa")
  on.exit(unlink(path), add = TRUE)
  writeLines(
    c(">seq1", "AAAA", ">seq2", "CCCC", ">seq3", "GGGG", ">seq4", "TTTT"),
    path
  )
  idx <- create_index(path, type = "fasta")

  out <- extract_sequences_dnastringset(
    idx,
    c(4L, 2L),
    mode = "sequential"
  )
  expect_equal(names(out), c("seq4", "seq2"))
  expect_equal(unname(as.character(out)), c("TTTT", "CCCC"))
})

test_that("extract_sequences_dnastringset streaming covers gzip FASTQ index", {
  skip_if_not_installed("Biostrings")
  path <- tempfile(fileext = ".fq.gz")
  on.exit(unlink(path), add = TRUE)
  make_fastq_gz(path)
  idx <- create_index(path, type = "fastq")

  out <- extract_sequences_dnastringset(
    idx,
    c(4L, 2L, 1L),
    mode = "sequential"
  )
  expect_equal(names(out), c("r4", "r2", "r1"))
  expect_equal(unname(as.character(out)), c("NANA", "TTAA", "ACGT"))
})

test_that("extract_sequences_dnastringset streaming works without index", {
  skip_if_not_installed("Biostrings")
  path <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(path), add = TRUE)
  make_fasta_gz(path)

  out <- extract_sequences_dnastringset(
    index = NULL,
    seq_idx = c(2L, 4L),
    file = path,
    type = "fasta"
  )
  expect_equal(names(out), c("seq2", "seq4"))
  expect_equal(unname(as.character(out)), c("CCCC", "TTTT"))
})

test_that("extract_sequences_dnastringset streaming honors renumber and empty input", {
  skip_if_not_installed("Biostrings")
  path <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(path), add = TRUE)
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")

  zero <- extract_sequences_dnastringset(
    idx,
    c(3L, 1L),
    mode = "sequential",
    renumber = "zero_based"
  )
  expect_equal(names(zero), c("0", "1"))

  one <- extract_sequences_dnastringset(
    idx,
    c(3L, 1L),
    mode = "sequential",
    renumber = "one_based"
  )
  expect_equal(names(one), c("1", "2"))

  empty <- extract_sequences_dnastringset(
    idx,
    integer(0L),
    mode = "sequential"
  )
  expect_s4_class(empty, "DNAStringSet")
  expect_identical(length(empty), 0L)
})

test_that("extract_sequences_dnastringset accepts list seq_idx", {
  skip_if_not_installed("Biostrings")
  path <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(path), add = TRUE)
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")

  parts <- list(c(1L, 2L), integer(), c(4L))
  out <- extract_sequences_dnastringset(idx, parts, mode = "sequential")
  expect_type(out, "list")
  expect_length(out, 3L)
  expect_equal(unname(as.character(out[[1L]])), c("AAAA", "CCCC"))
  expect_identical(length(out[[2L]]), 0L)
  expect_equal(unname(as.character(out[[3L]])), "TTTT")
})

test_that("extract_sequences_dnastringset validates chunk_chars", {
  skip_if_not_installed("Biostrings")
  path <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(path), add = TRUE)
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")

  expect_error(
    extract_sequences_dnastringset(idx, 1L, chunk_chars = -1),
    "`chunk_chars` must be"
  )
  expect_error(
    extract_sequences_dnastringset(idx, 1L, chunk_chars = NA_real_),
    "`chunk_chars` must be"
  )
  expect_error(
    extract_sequences_dnastringset(idx, 1L, chunk_chars = Inf),
    "`chunk_chars` must be"
  )
  expect_error(
    extract_sequences_dnastringset(idx, 1L, chunk_chars = c(1, 2)),
    "`chunk_chars` must be"
  )
})

test_that("extract_sequences_dnastringset rejects file override length mismatch", {
  skip_if_not_installed("Biostrings")
  p1 <- tempfile(fileext = ".fa.gz")
  p2 <- tempfile(fileext = ".fa.gz")
  alt <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(c(p1, p2, alt)), add = TRUE)
  make_fasta_gz(p1)
  make_fasta_gz(p2)
  make_fasta_gz(alt)
  idx <- create_index(c(p1, p2), type = "fasta")

  expect_error(
    extract_sequences_dnastringset(idx, 1L, file = alt),
    "same number of files"
  )
})

test_that("extract_sequences_dnastringset with mode='indexed' requires non-NULL index", {
  skip_if_not_installed("Biostrings")
  path <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(path), add = TRUE)
  make_fasta_gz(path)

  expect_error(
    extract_sequences_dnastringset(
      index = NULL,
      seq_idx = 1L,
      file = path,
      mode = "indexed"
    ),
    "requires a non-NULL `index`"
  )
})
