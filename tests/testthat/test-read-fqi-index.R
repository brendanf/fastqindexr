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
