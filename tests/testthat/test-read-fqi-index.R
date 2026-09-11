fixture_path <- function(...) {
  testthat::test_path("fixtures", ...)
}

test_that("read_fqi_index reads fixture with explicit files", {
  fqi <- fixture_path("cli_fixture.fastq.gz.fqi")
  fq <- fixture_path("cli_fixture.fastq.gz")

  idx <- read_fqi_index(fqi_path = fqi, files = fq, type = "fastq")

  expect_s3_class(idx, "fastqindexr_index")
  expect_s3_class(idx, "fastqindexr_gzip_index")
  expect_identical(idx$format, "fastq")
  expect_equal(idx$n_records, 4)
  expect_equal(unname(idx$file_record_offsets), c(0, 4))
})

test_that("read_fqi_index deduces files and type", {
  fqi <- fixture_path("cli_fixture.fastq.gz.fqi")
  idx <- read_fqi_index(fqi_path = fqi, files = NULL, type = "auto")

  expect_identical(idx$format, "fastq")
  expect_equal(idx$n_records, 4)
})

test_that("read_fqi_index warns when deduced source file is missing", {
  src_fqi <- fixture_path("cli_fixture.fastq.gz.fqi")
  tmp_fqi <- tempfile(fileext = ".fastq.gz.fqi")
  file.copy(src_fqi, tmp_fqi, overwrite = TRUE)

  expect_warning(
    idx <- read_fqi_index(fqi_path = tmp_fqi, files = NULL, type = "fastq"),
    "Deduced source file"
  )
  expect_s3_class(idx, "fastqindexr_index")
  expect_s3_class(idx, "fastqindexr_gzip_index")
})

test_that("read_fqi_index warns when auto type cannot be inferred", {
  # One `read_fqi_index` call emits three warnings. testthat 3e links
  # them on a parent chain, so a pipe of `expect_warning` calls can
  # assert each; any other warning is not matched and is reported.
  src_fqi <- fixture_path("cli_fixture.fastq.gz.fqi")
  tmp_fqi <- tempfile(fileext = ".fastq.gz.fqi")
  file.copy(src_fqi, tmp_fqi, overwrite = TRUE)

  (idx <- read_fqi_index(fqi_path = tmp_fqi, files = NULL, type = "auto")) |>
    expect_warning("Deduced source file") |>
    expect_warning("Could not determine type from file") |>
    expect_warning("defaulting")
  expect_s3_class(idx, "fastqindexr_index")
  expect_s3_class(idx, "fastqindexr_gzip_index")
  expect_identical(idx$format, "fastq")
})

test_that("extract_sequences accepts .fqi path with explicit file", {
  fqi <- fixture_path("cli_fixture.fastq.gz.fqi")
  fq <- fixture_path("cli_fixture.fastq.gz")
  idx <- create_index(fq, type = "fastq")

  from_obj <- extract_sequences(idx, c(4, 2, 1))
  from_fqi <- extract_sequences(fqi, c(4, 2, 1), file = fq)

  expect_identical(from_fqi, from_obj)
})

test_that("extract_sequences accepts .fqi path with deduced file", {
  fqi <- fixture_path("cli_fixture.fastq.gz.fqi")
  out <- extract_sequences(fqi, c(3, 1), file = NULL)

  expect_equal(out$seq_id, c("r3", "r1"))
  expect_equal(out$seq, c("GGGG", "ACGT"))
})

test_that("read_fqi_index supports vector fqi_path and aligned files", {
  fqi <- c(
    fixture_path("cli_fixture.fastq.gz.fqi"),
    fixture_path("cli_fixture_part2.fastq.gz.fqi")
  )
  fq <- c(
    fixture_path("cli_fixture.fastq.gz"),
    fixture_path("cli_fixture_part2.fastq.gz")
  )

  idx <- read_fqi_index(fqi_path = fqi, files = fq, type = "fastq")
  expect_equal(idx$n_records, 6)
  expect_equal(unname(idx$file_record_offsets), c(0, 4, 6))

  out <- extract_sequences(idx, c(5, 2, 6))
  expect_equal(out$seq_id, c("r5", "r2", "r6"))
})

test_that("read_fqi_index errors when files length disagrees with fqi_path", {
  fqi <- c(
    fixture_path("cli_fixture.fastq.gz.fqi"),
    fixture_path("cli_fixture_part2.fastq.gz.fqi")
  )
  fq <- fixture_path("cli_fixture.fastq.gz")
  expect_error(
    read_fqi_index(fqi_path = fqi, files = fq, type = "fastq"),
    "same length as `fqi_path`"
  )
})

test_that("detect_type_from_gz returns NA for empty/blank/unknown content", {
  detect <- asNamespace("fastqindexr")$detect_type_from_gz

  empty_path <- tempfile(fileext = ".gz")
  on.exit(unlink(empty_path), add = TRUE)
  con <- gzfile(empty_path, "wb")
  close(con)
  expect_true(is.na(detect(empty_path)))

  blank_path <- tempfile(fileext = ".gz")
  on.exit(unlink(blank_path), add = TRUE)
  blank_con <- gzfile(blank_path, "wt")
  writeLines(c("", "", ""), blank_con)
  close(blank_con)
  expect_true(is.na(detect(blank_path)))

  unknown_path <- tempfile(fileext = ".gz")
  on.exit(unlink(unknown_path), add = TRUE)
  unk_con <- gzfile(unknown_path, "wt")
  writeLines(c("not a header", "ACGT"), unk_con)
  close(unk_con)
  expect_true(is.na(detect(unknown_path)))

  fasta_path <- tempfile(fileext = ".gz")
  on.exit(unlink(fasta_path), add = TRUE)
  fa_con <- gzfile(fasta_path, "wt")
  writeLines(c(">s1", "AAAA"), fa_con)
  close(fa_con)
  expect_identical(detect(fasta_path), "fasta")
})

test_that("read_fqi_index errors when type='auto' detects mixed FASTA/FASTQ", {
  resolve <- asNamespace("fastqindexr")$resolve_type_for_files

  fa_gz <- tempfile(fileext = ".fa.gz")
  fq_gz <- tempfile(fileext = ".fq.gz")
  on.exit(unlink(c(fa_gz, fq_gz)), add = TRUE)
  fa_con <- gzfile(fa_gz, "wt")
  writeLines(c(">s1", "AAAA"), fa_con)
  close(fa_con)
  fq_con <- gzfile(fq_gz, "wt")
  writeLines(c("@r1", "AAAA", "+", "!!!!"), fq_con)
  close(fq_con)

  expect_error(
    resolve("auto", c(fa_gz, fq_gz)),
    "mixed FASTA/FASTQ"
  )
})

test_that("resolve_type_for_files returns explicit type without scanning", {
  resolve <- asNamespace("fastqindexr")$resolve_type_for_files
  expect_identical(resolve("fasta", c("nope1", "nope2")), "fasta")
  expect_identical(resolve("fastq", c("nope1")), "fastq")
})
