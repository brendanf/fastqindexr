# validate_input_files ---------------------------------------------------------

test_that("create_index errors when files vector is empty", {
  expect_error(create_index(character(), type = "fasta"), "at least one path")
})

test_that("create_index errors when a path does not exist", {
  missing <- file.path(tempdir(), "definitely-not-a-real-file.fa.gz")
  expect_error(create_index(missing, type = "fasta"), "Missing file")
})

test_that("create_index validates index_stride_bytes", {
  path <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(path), add = TRUE)
  make_fasta_gz(path)

  expect_error(
    create_index(path, type = "fasta", index_stride_bytes = "x"),
    "`index_stride_bytes` must be"
  )
  expect_error(
    create_index(path, type = "fasta", index_stride_bytes = NA_real_),
    "`index_stride_bytes` must be"
  )
  expect_error(
    create_index(path, type = "fasta", index_stride_bytes = c(1, 2)),
    "`index_stride_bytes` must be"
  )
  # Numeric scalar is accepted (currently ignored).
  idx <- create_index(path, type = "fasta", index_stride_bytes = 100)
  expect_s3_class(idx, "fastqindexr_index")
})

# assign_index_classes ---------------------------------------------------------

test_that("assign_index_classes rejects mixed compression in payload", {
  internal <- asNamespace("fastqindexr")$assign_index_classes
  bad <- list(
    files = c("a", "b"),
    format = "fasta",
    n_records = 0,
    file_record_offsets = c(0, 0, 0),
    index_payload = list(),
    file_compression = c("gzip", "plain")
  )
  expect_error(internal(bad), "mixed file compression")
})

test_that("assign_index_classes errors when subclass cannot be inferred", {
  internal <- asNamespace("fastqindexr")$assign_index_classes
  no_compression <- list(
    files = "a",
    format = "fasta",
    n_records = 0,
    file_record_offsets = c(0, 0),
    index_payload = list(),
    file_compression = character()
  )
  expect_error(
    internal(no_compression),
    "Could not determine index subclass"
  )

  unknown_compression <- list(
    files = "a",
    format = "fasta",
    n_records = 0,
    file_record_offsets = c(0, 0),
    index_payload = list(),
    file_compression = "unknown"
  )
  expect_error(
    internal(unknown_compression),
    "Could not determine index subclass"
  )
})

# validate_index ---------------------------------------------------------------

test_that("extract_sequences rejects non-fastqindexr_index objects", {
  expect_error(
    extract_sequences(structure(list(), class = "not_an_index"), 1),
    "must be NULL, a fastqindexr_index"
  )
})

test_that("validate_index detects missing required fields", {
  internal <- asNamespace("fastqindexr")$validate_index
  bad <- structure(
    list(format = "fasta"),
    class = c("fastqindexr_gzip_index", "fastqindexr_index", "list")
  )
  expect_error(internal(bad), "Missing: ")
})

test_that("validate_index rejects malformed file_compression", {
  internal <- asNamespace("fastqindexr")$validate_index
  base <- list(
    files = "a",
    format = "fasta",
    n_records = 0,
    file_record_offsets = c(0, 0),
    index_payload = list()
  )

  bad_no_comp <- structure(
    c(base, list(file_compression = character())),
    class = c("fastqindexr_gzip_index", "fastqindexr_index", "list")
  )
  expect_error(internal(bad_no_comp), "homogeneous")

  bad_mixed <- structure(
    c(base, list(file_compression = c("gzip", "plain"))),
    class = c("fastqindexr_gzip_index", "fastqindexr_index", "list")
  )
  expect_error(internal(bad_mixed), "homogeneous")

  bad_unknown <- structure(
    c(base, list(file_compression = "weird")),
    class = c("fastqindexr_gzip_index", "fastqindexr_index", "list")
  )
  expect_error(internal(bad_unknown), "homogeneous")
})

test_that("validate_index enforces gzip subclass <-> compression coherence", {
  internal <- asNamespace("fastqindexr")$validate_index
  base <- list(
    files = "a",
    format = "fasta",
    n_records = 0,
    file_record_offsets = c(0, 0),
    index_payload = list()
  )

  gzip_with_plain <- structure(
    c(base, list(file_compression = "plain")),
    class = c("fastqindexr_gzip_index", "fastqindexr_index", "list")
  )
  expect_error(
    internal(gzip_with_plain),
    "gzip subclass requires gzip compression"
  )

  plain_with_gzip <- structure(
    c(base, list(file_compression = "gzip")),
    class = c("fastqindexr_plain_index", "fastqindexr_index", "list")
  )
  expect_error(
    internal(plain_with_gzip),
    "plain subclass requires plain compression"
  )
})

# resolve_extract_index --------------------------------------------------------

test_that("extract_sequences rejects index that is neither object nor path", {
  expect_error(
    extract_sequences(index = 42, seq_idx = 1, file = NULL),
    "must be NULL, a fastqindexr_index"
  )
})

# validate_seq_idx -------------------------------------------------------------

test_that("extract_sequences rejects non-numeric seq_idx", {
  path <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(path), add = TRUE)
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")
  expect_error(
    extract_sequences(idx, seq_idx = c("a", "b")),
    "must be numeric/integer-like"
  )
})

# extract_sequences mode/file overrides ----------------------------------------

test_that("extract_sequences with mode='indexed' requires non-NULL index", {
  path <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(path), add = TRUE)
  make_fasta_gz(path)
  expect_error(
    extract_sequences(
      index = NULL,
      seq_idx = 1,
      file = path,
      mode = "indexed"
    ),
    "requires a non-NULL `index`"
  )
})

test_that("extract_sequences rejects file override of mismatched length", {
  p1 <- tempfile(fileext = ".fa.gz")
  p2 <- tempfile(fileext = ".fa.gz")
  alt <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(c(p1, p2, alt)), add = TRUE)
  make_fasta_gz(p1)
  make_fasta_gz(p2)
  make_fasta_gz(alt)
  idx <- create_index(c(p1, p2), type = "fasta")
  expect_error(
    extract_sequences(idx, c(1L, 2L), file = alt),
    "same number of files"
  )
})

# extract_sequences_to_file argument validation --------------------------------

test_that("extract_sequences_to_file validates append, outfile, collapse", {
  path <- tempfile(fileext = ".fa.gz")
  out <- tempfile(fileext = ".fa")
  on.exit(unlink(c(path, out)), add = TRUE)
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")

  expect_error(
    extract_sequences_to_file(idx, 1, outfile = out, append = NA),
    "`append` must be"
  )
  expect_error(
    extract_sequences_to_file(idx, 1, outfile = out, append = "yes"),
    "`append` must be"
  )
  expect_error(
    extract_sequences_to_file(idx, 1, outfile = character()),
    "`outfile` must be"
  )
  expect_error(
    extract_sequences_to_file(idx, 1, outfile = ""),
    "`outfile` must be"
  )
  expect_error(
    extract_sequences_to_file(
      idx,
      1,
      outfile = out,
      collapse_sequence_lines = NA
    ),
    "`collapse_sequence_lines` must be"
  )
  expect_error(
    extract_sequences_to_file(
      idx,
      1,
      outfile = out,
      compress = NA
    ),
    "`compress` must be"
  )
})

test_that("extract_sequences_to_file with mode='indexed' requires non-NULL index", {
  path <- tempfile(fileext = ".fa.gz")
  out <- tempfile(fileext = ".fa")
  on.exit(unlink(c(path, out)), add = TRUE)
  make_fasta_gz(path)

  expect_error(
    extract_sequences_to_file(
      index = NULL,
      seq_idx = 1,
      file = path,
      outfile = out,
      mode = "indexed"
    ),
    "requires a non-NULL `index`"
  )
})

test_that("extract_sequences_to_file rejects file override of mismatched length", {
  p1 <- tempfile(fileext = ".fa.gz")
  p2 <- tempfile(fileext = ".fa.gz")
  alt <- tempfile(fileext = ".fa.gz")
  out <- tempfile(fileext = ".fa")
  on.exit(unlink(c(p1, p2, alt, out)), add = TRUE)
  make_fasta_gz(p1)
  make_fasta_gz(p2)
  make_fasta_gz(alt)
  idx <- create_index(c(p1, p2), type = "fasta")
  expect_error(
    extract_sequences_to_file(idx, 1, file = alt, outfile = out),
    "same number of files"
  )
})
