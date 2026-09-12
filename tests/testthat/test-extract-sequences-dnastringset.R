test_that("extract_sequences_dnastringset preserves order and duplicates", {
  path <- tempfile(fileext = ".fa.gz")
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")

  out <- extract_sequences_dnastringset(idx, c(3, 1, 3))
  expect_s4_class(out, "DNAStringSet")
  expect_equal(names(out), c("seq3", "seq1", "seq3"))
  expect_equal(unname(as.character(out)), c("GGGG", "AAAA", "GGGG"))
})

test_that("extract_sequences_dnastringset renumber modes are deterministic", {
  path <- tempfile(fileext = ".fa.gz")
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")
  ids <- c(4, 2, 1)

  zero <- extract_sequences_dnastringset(idx, ids, renumber = "zero_based")
  one <- extract_sequences_dnastringset(idx, ids, renumber = "one_based")
  expect_equal(names(zero), c("0", "1", "2"))
  expect_equal(names(one), c("1", "2", "3"))
  expect_equal(unname(as.character(zero)), unname(as.character(one)))
})

test_that("extract_sequences_dnastringset concatenates multiple chunks", {
  path <- tempfile(fileext = ".fa.gz")
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")
  # 4-char sequences; a tiny limit forces more than one DNAStringSet chunk.
  ids <- c(3, 1, 3, 2, 4)
  out <- extract_sequences_dnastringset(
    idx,
    ids,
    chunk_chars = 6,
    renumber = "none"
  )
  expect_equal(names(out), c("seq3", "seq1", "seq3", "seq2", "seq4"))
  expect_equal(
    unname(as.character(out)),
    c("GGGG", "AAAA", "GGGG", "CCCC", "TTTT")
  )
})

test_that("extract_sequences_dnastringset works for FASTQ index and empty input", {
  path <- tempfile(fileext = ".fq.gz")
  make_fastq_gz(path)
  idx <- create_index(path, type = "fastq")

  out <- extract_sequences_dnastringset(idx, c(4, 2))
  expect_equal(names(out), c("r4", "r2"))
  expect_equal(unname(as.character(out)), c("NANA", "TTAA"))

  empty <- extract_sequences_dnastringset(idx, integer(0))
  expect_s4_class(empty, "DNAStringSet")
  expect_identical(length(empty), 0L)
})
