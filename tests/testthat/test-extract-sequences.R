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

test_that("extract_sequences restores native pointer from serialized payload", {
  path <- tempfile(fileext = ".fa.gz")
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")

  idx$._cache$index_ptr <- NULL
  out <- extract_sequences(idx, c(4, 2))
  expect_equal(out$seq_id, c("seq4", "seq2"))
})

test_that("serialized index object works after readRDS", {
  path <- tempfile(fileext = ".fa.gz")
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")
  rds <- tempfile(fileext = ".rds")
  saveRDS(idx, rds)

  restored <- readRDS(rds)
  out <- extract_sequences(restored, c(3, 1))
  expect_equal(out$seq_id, c("seq3", "seq1"))
})

test_that("extract_sequences return = list for FASTA and FASTQ", {
  fa_path <- tempfile(fileext = ".fa.gz")
  make_fasta_gz(fa_path)
  fa_idx <- create_index(fa_path, type = "fasta")
  l_fa <- extract_sequences(fa_idx, c(3, 1, 2), return = "list")
  expect_type(l_fa, "list")
  expect_setequal(names(l_fa), c("seq_id", "seq"))
  expect_equal(l_fa$seq_id, c("seq3", "seq1", "seq2"))
  expect_equal(l_fa$seq, c("GGGG", "AAAA", "CCCC"))

  fq_path <- tempfile(fileext = ".fq.gz")
  make_fastq_gz(fq_path)
  fq_idx <- create_index(fq_path, type = "fastq")
  l_fq <- extract_sequences(fq_idx, c(4, 2), return = "list")
  expect_setequal(names(l_fq), c("seq_id", "seq", "qual"))
  expect_equal(l_fq$seq_id, c("r4", "r2"))
  expect_equal(l_fq$qual, c("%%%%", "####"))
})

test_that("extract_sequences return = seq with duplicate names and order", {
  path <- tempfile(fileext = ".fa.gz")
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")
  s <- extract_sequences(idx, c(3, 1, 3), return = "seq")
  expect_type(s, "character")
  expect_equal(unname(s), c("GGGG", "AAAA", "GGGG"))
  expect_equal(names(s), c("seq3", "seq1", "seq3"))
})

test_that("extract_sequences return = seq for FASTQ matches data.frame seq column", {
  path <- tempfile(fileext = ".fq.gz")
  make_fastq_gz(path)
  idx <- create_index(path, type = "fastq")
  ids <- c(4, 2, 1)
  as_df <- extract_sequences(idx, ids, return = "data.frame")
  as_seq <- extract_sequences(idx, ids, return = "seq")
  expect_equal(unname(as_seq), as_df$seq)
  expect_equal(names(as_seq), as_df$seq_id)
})

test_that("extract_sequences empty extract respects return mode", {
  path <- tempfile(fileext = ".fa.gz")
  make_fasta_gz(path)
  fa <- create_index(path, type = "fasta")
  expect_equal(
    extract_sequences(fa, integer(0L), return = "data.frame"),
    data.frame(
      seq_id = character(),
      seq = character(),
      stringsAsFactors = FALSE
    )
  )
  el <- extract_sequences(fa, integer(0L), return = "list")
  expect_type(el, "list")
  expect_setequal(names(el), c("seq_id", "seq"))
  expect_equal(
    extract_sequences(fa, integer(0L), return = "seq"),
    setNames(character(0L), character(0L))
  )

  fq_path <- tempfile(fileext = ".fq.gz")
  make_fastq_gz(fq_path)
  fq <- create_index(fq_path, type = "fastq")
  lq <- extract_sequences(fq, integer(0L), return = "list")
  expect_setequal(names(lq), c("seq_id", "seq", "qual"))
})

test_that("extraction matches for sorted unique vs same multiset permuted", {
  tmp <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(tmp), add = TRUE)
  make_benchmark_fasta(tmp, n = 500L, width = 40L)
  idx <- create_index(tmp, type = "fasta")
  k <- 200L
  set.seed(42L)
  ids <- sample.int(500L, k, replace = FALSE)
  a <- extract_sequences(idx, ids)
  b <- extract_sequences(idx, sort(ids))
  ord <- order(ids)
  expect_equal(b$seq_id, a$seq_id[ord])
  expect_equal(b$seq, a$seq[ord])
})

test_that("extractors honor region merge tuning options", {
  tmp <- tempfile(fileext = ".fa.gz")
  out_file <- tempfile(fileext = ".fa")
  on.exit(unlink(c(tmp, out_file)), add = TRUE)
  make_benchmark_fasta(tmp, n = 200L, width = 30L)
  idx <- create_index(tmp, type = "fasta")
  ids <- c(190L, 2L, 155L, 99L, 155L, 11L)

  baseline <- extract_sequences(idx, ids, return = "list")
  old_opts <- options(
    fastqindexr.max_bridge_gap = 0,
    fastqindexr.max_region_bytes = 8
  )
  on.exit(options(old_opts), add = TRUE)

  tuned <- extract_sequences(idx, ids, return = "list")
  expect_equal(tuned$seq_id, baseline$seq_id)
  expect_equal(tuned$seq, baseline$seq)

  extract_sequences_to_file(idx, ids, outfile = out_file, append = FALSE)
  lines <- readLines(out_file, warn = FALSE)
  expect_equal(
    lines[seq(1L, length(lines), by = 2L)],
    paste0(">", baseline$seq_id)
  )

  dna <- extract_sequences_dnastringset(idx, ids)
  expect_equal(names(dna), baseline$seq_id)
  expect_equal(unname(as.character(dna)), baseline$seq)
})

test_that("region merge tuning options validate cleanly", {
  path <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(path), add = TRUE)
  make_fasta_gz(path)
  idx <- create_index(path, type = "fasta")

  old_opts <- options(
    fastqindexr.max_bridge_gap = NA_real_,
    fastqindexr.max_region_bytes = 2147483647,
    fastqindexr.extract_mode = "indexed",
    fastqindexr.extract_diagnostics = FALSE
  )
  on.exit(options(old_opts), add = TRUE)
  expect_error(
    extract_sequences(idx, 1L),
    "fastqindexr.max_bridge_gap"
  )

  options(fastqindexr.max_bridge_gap = 1L, fastqindexr.max_region_bytes = Inf)
  expect_error(
    extract_sequences(idx, 1L),
    "fastqindexr.max_region_bytes"
  )

  options(
    fastqindexr.max_bridge_gap = 1L,
    fastqindexr.max_region_bytes = 1000L,
    fastqindexr.extract_mode = "bad"
  )
  expect_error(
    extract_sequences(idx, 1L),
    "fastqindexr.extract_mode"
  )

  options(
    fastqindexr.max_bridge_gap = 1L,
    fastqindexr.max_region_bytes = 1000L,
    fastqindexr.extract_mode = "indexed",
    fastqindexr.extract_diagnostics = 1
  )
  expect_error(
    extract_sequences(idx, 1L),
    "fastqindexr.extract_diagnostics"
  )
})

test_that("sequential_only mode preserves extraction outputs", {
  tmp <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(tmp), add = TRUE)
  make_benchmark_fasta(tmp, n = 500L, width = 20L)
  idx <- create_index(tmp, type = "fasta")
  set.seed(99L)
  ids <- sample.int(500L, 120L, replace = TRUE)

  old_opts <- options(
    fastqindexr.max_bridge_gap = 64,
    fastqindexr.max_region_bytes = 2147483647,
    fastqindexr.extract_mode = "indexed",
    fastqindexr.extract_diagnostics = FALSE
  )
  on.exit(options(old_opts), add = TRUE)

  baseline <- extract_sequences(idx, ids, return = "list")
  options(fastqindexr.extract_mode = "sequential_only")
  got <- extract_sequences(idx, ids, return = "list")
  expect_equal(got, baseline)
})

test_that("diagnostics attributes are exposed when enabled", {
  tmp <- tempfile(fileext = ".fa.gz")
  out <- tempfile(fileext = ".fa")
  on.exit(unlink(c(tmp, out)), add = TRUE)
  make_benchmark_fasta(tmp, n = 200L, width = 30L)
  idx <- create_index(tmp, type = "fasta")

  old_opts <- options(
    fastqindexr.max_bridge_gap = 64,
    fastqindexr.max_region_bytes = 2147483647,
    fastqindexr.extract_mode = "indexed",
    fastqindexr.extract_diagnostics = TRUE
  )
  on.exit(options(old_opts), add = TRUE)

  res <- extract_sequences(idx, c(1L, 50L, 150L), return = "list")
  diag <- attr(res, "fastqindexr_diagnostics", exact = TRUE)
  expect_true(is.list(diag))
  expect_setequal(
    names(diag),
    c(
      "regions_planned",
      "extract_attempts",
      "extract_failures",
      "fallback_invocations",
      "fallback_records"
    )
  )

  out_path <- extract_sequences_to_file(
    idx,
    seq_idx = c(1L, 50L, 150L),
    outfile = out
  )
  out_diag <- attr(out_path, "fastqindexr_diagnostics", exact = TRUE)
  expect_true(is.list(out_diag))
  expect_true(out_diag$extract_attempts >= 0)
})

test_that("extract_sequences output is stable across density and region-merge grid", {
  tmp <- tempfile(fileext = ".fa.gz")
  on.exit(unlink(tmp), add = TRUE)
  make_benchmark_fasta(tmp, n = 800L, width = 25L)
  idx <- create_index(tmp, type = "fasta")

  dense <- 200:219L
  medium <- c(200:205, 260:265, 320:325)
  sparse <- c(5L, 400L, 790L)
  patterns <- list(
    dense = dense,
    medium = medium,
    sparse = sparse
  )

  opt_grid <- expand.grid(
    max_bridge_gap = c(0L, 64L, 256L),
    max_region_bytes = c(50000L, 500000L, 2147483647L),
    stringsAsFactors = FALSE
  )

  old_opts <- options(
    fastqindexr.max_bridge_gap = 64,
    fastqindexr.max_region_bytes = 2147483647,
    fastqindexr.extract_mode = "indexed",
    fastqindexr.extract_diagnostics = FALSE
  )
  on.exit(options(old_opts), add = TRUE)

  for (nm in names(patterns)) {
    ids <- patterns[[nm]]
    options(
      fastqindexr.max_bridge_gap = 64,
      fastqindexr.max_region_bytes = 2147483647,
      fastqindexr.extract_mode = "indexed",
      fastqindexr.extract_diagnostics = FALSE
    )
    baseline <- extract_sequences(idx, ids, return = "list")
    for (k in seq_len(nrow(opt_grid))) {
      options(
        fastqindexr.max_bridge_gap = opt_grid$max_bridge_gap[k],
        fastqindexr.max_region_bytes = opt_grid$max_region_bytes[k],
        fastqindexr.extract_mode = "indexed",
        fastqindexr.extract_diagnostics = FALSE
      )
      got <- extract_sequences(idx, ids, return = "list")
      expect_equal(got, baseline, info = sprintf("%s grid row %d", nm, k))
    }
  }
})
