#' @noRd
validate_index <- function(index) {
  if (!inherits(index, "fastqindexr_index")) {
    stop("`index` must be a fastqindexr_index object.", call. = FALSE)
  }
  required <- c(
    "files",
    "format",
    "n_records",
    "file_record_offsets",
    "index_payload"
  )
  missing <- setdiff(required, names(index))
  if (length(missing) > 0L) {
    missing_fields <- paste(missing, collapse = ", ")
    stop(
      sprintf("Malformed index object. Missing: %s", missing_fields),
      call. = FALSE
    )
  }
  index
}

#' @noRd
resolve_extract_index <- function(index, file) {
  if (inherits(index, "fastqindexr_index")) {
    return(list(index = index, file = file))
  }
  if (is.character(index)) {
    loaded <- read_fqi_index(fqi_path = index, files = file, type = "auto")
    return(list(index = loaded, file = NULL))
  }
  stop(
    "`index` must be a fastqindexr_index object or path(s) to `.fqi` file(s).",
    call. = FALSE
  )
}

#' @noRd
# nolint start: object_usage_linter
ensure_live_index_ptr <- function(index) {
  cache <- index$._cache
  if (!is.environment(cache)) {
    cache <- new.env(parent = emptyenv())
    index$._cache <- cache
  }

  ptr <- cache$index_ptr
  if (!isTRUE(cpp_index_ptr_is_valid(ptr))) {
    ptr <- cpp_restore_index_ptr(index$index_payload)
    cache$index_ptr <- ptr
  }

  list(index = index, ptr = ptr)
}
# nolint end

#' Extract sequence records by ID from indexed gzipped FASTA or FASTQ
#'
#' Returns rows in **the same order as `seq_idx`**, including duplicate IDs. IDs
#' are
#' **1-based** record indices in the logical concatenated stream defined when
#' the index was built (see [create_index()]).
#'
#' @param index Either a `fastqindexr_index` object (from [create_index()] or
#'   [read_fqi_index()]) or one or more `.fqi` file paths.
#' @param seq_idx Numeric vector of record IDs (positive whole numbers). Values
#'   are coerced via [as.numeric()]; `NA` and values outside `1` &hellip; `n_records`
#'   are errors.
#' @param file Optional character vector of file paths overriding those
#'   stored in `index$files` for object input.
#'   For `.fqi` path input, this provides indexed gz file path(s) passed to
#'   [read_fqi_index()]. If omitted for `.fqi` input, file paths are deduced from
#'   the `.fqi` names.
#' @param return In-memory return shape: `"data.frame"` (the default), `"list"`, or
#'   `"seq"`. With `"data.frame"`, the result is as before. With `"list"`, returns
#'   a list with `seq_id` and `seq` (and `qual` for FASTQ). With `"seq"`, returns
#'   a character vector of sequences with `seq_id` as names, preserving order
#'   and allowing duplicate names. For indexed FASTQ, `return = "seq"` skips
#'   reading quality lines from the source.
#'
#' @return
#' When `return = "data.frame"` (the default), a `data.frame`:
#' \describe{
#'   \item{FASTA}{Columns `seq_id` and `seq`.}
#'   \item{FASTQ}{Columns `seq_id`, `seq`, and `qual`.}
#' }
#' For `length(seq_idx) == 0`, an empty `data.frame` with the appropriate columns
#' is returned.
#'
#' When `return = "list"`, a `list` with `seq_id`, `seq`, and for FASTQ also
#' `qual` (omitted for FASTA).
#'
#' When `return = "seq"`, a character vector of `seq` only, with names
#' set from `seq_id` (duplicates allowed).
#' Empty extract requests use empty containers of the corresponding type.
#'
#' @section FASTA limitation:
#' Each record must consist of a header line plus **one** sequence line.
#'
#' @seealso [create_index()]
#'
#' @examples
#' path <- tempfile(fileext = ".fastq.gz")
#' con <- gzfile(path, "wt")
#' writeLines(c("@r1", "ACGT", "+", "!!!!", "@r2", "TTAA", "+", "####"), con)
#' close(con)
#' idx <- create_index(path, type = "fastq")
#' extract_sequences(idx, seq_idx = c(2, 2, 1))
#' unlink(path)
#'
#' @export
extract_sequences <- function(
  index,
  seq_idx,
  file = NULL,
  return = c("data.frame", "list", "seq")
) {
  return <- match.arg(return)
  resolved <- resolve_extract_index(index, file)
  index <- resolved$index
  file <- resolved$file
  index <- validate_index(index)
  live_index <- ensure_live_index_ptr(index)
  index <- live_index$index
  index_ptr <- live_index$ptr

  if (length(seq_idx) < 1L) {
    if (return == "list") {
      if (identical(index$format, "fastq")) {
        return(
          list(
            seq_id = character(),
            seq = character(),
            qual = character()
          )
        )
      }
      return(list(seq_id = character(), seq = character()))
    }
    if (return == "seq") {
      return(structure(character(0L), names = character(0L)))
    }
    if (identical(index$format, "fastq")) {
      return(data.frame(
        seq_id = character(),
        seq = character(),
        qual = character(),
        stringsAsFactors = FALSE
      ))
    }
    return(data.frame(
      seq_id = character(),
      seq = character(),
      stringsAsFactors = FALSE
    ))
  }

  if (!is.numeric(seq_idx)) {
    stop("`seq_idx` must be numeric/integer-like.", call. = FALSE)
  }
  if (any(is.na(seq_idx) | seq_idx < 1 | seq_idx != floor(seq_idx))) {
    stop(
      "`seq_idx` must contain positive whole numbers (1-based).",
       call. = FALSE
     )
  }

  n_records <- as.numeric(index$n_records)
  if (any(seq_idx > n_records)) {
    stop(
      sprintf("Some seq_idx exceed available records (%s).", n_records),
      call. = FALSE
    )
  }

  files <- if (is.null(file)) {
    as.character(index$files)
  } else {
    files <- validate_input_files(file)
    if (length(files) != length(index$files)) {
      stop(
        paste0(
          "Override `file` must contain the same number of files as ",
          "index$files."
        ),
        call. = FALSE
      )
    }
    files
  }

  is_fastq <- identical(index$format, "fastq")
  include_qual <- is_fastq && (return %in% c("data.frame", "list"))

  # nolint nextline: object_usage_linter
  result <- cpp_extract_sequences(
    files = files,
    type = index$format,
    ids_zero_based = as.numeric(seq_idx - 1),
    index_ptr_sexp = index_ptr,
    include_qual = include_qual
  )

  if (return == "data.frame") {
    if (is_fastq) {
      return(
        data.frame(
          seq_id = result$seq_id,
          seq = result$seq,
          qual = result$qual,
          stringsAsFactors = FALSE
        )
      )
    }
    return(
      data.frame(
        seq_id = result$seq_id,
        seq = result$seq,
        stringsAsFactors = FALSE
      )
    )
  }
  if (return == "list") {
    if (is_fastq) {
      return(
        list(
          seq_id = result$seq_id,
          seq = result$seq,
          qual = result$qual
        )
      )
    }
    return(list(seq_id = result$seq_id, seq = result$seq))
  }
  stats::setNames(result$seq, result$seq_id)
}

#' Extract sequence records by ID and stream directly to file
#'
#' Writes records in **the same order as `seq_idx`**, including duplicate IDs.
#' IDs are **1-based** record indices in the logical concatenated stream defined
#' when the index was built (see [create_index()]).
#'
#' @param index Either a `fastqindexr_index` object (from [create_index()] or
#'   [read_fqi_index()]) or one or more `.fqi` file paths.
#' @param seq_idx Numeric vector of record IDs (positive whole numbers). Values
#'   are coerced via [as.numeric()]; `NA` and values outside `1` ... `n_records`
#'   are errors.
#' @param file Optional character vector of file paths overriding those stored
#'   in `index$files` for object input. For `.fqi` path input, this provides
#'   indexed gz file path(s) passed to [read_fqi_index()]. If omitted for `.fqi`
#'   input, file paths are deduced from the `.fqi` names.
#' @param outfile Output file path.
#' @param type Output format: `"auto"` (default), `"fasta"`, or `"fastq"`.
#'   `"auto"` uses the indexed input format. FASTQ input can be emitted as FASTA;
#'   FASTA input cannot be emitted as FASTQ. When emitting FASTA from a FASTQ
#'   index, quality lines are not read from the source.
#' @param append Logical; append to existing `outfile` if `TRUE`, otherwise
#'   overwrite. When `append` is `FALSE` and `seq_idx` is **strictly
#'   increasing** (no duplicate IDs), extraction writes each record as it is
#'   read and does not build a full in-memory map of all unique records. When
#'   `append` is `TRUE` or the request is not strictly increasing, the
#'   implementation may buffer all unique records like [extract_sequences()].
#' @param compress Logical; if `TRUE`, write gzip-compressed output directly via
#'   zlib. Defaults to `endsWith(tolower(outfile), ".gz")`.
#'
#' @return Invisibly returns `outfile`.
#'
#' @section FASTA limitation:
#' Each record must consist of a header line plus **one** sequence line.
#'
#' @seealso [extract_sequences()], [create_index()]
#'
#' @export
extract_sequences_to_file <- function(
  index,
  seq_idx,
  file = NULL,
  outfile,
  type = c("auto", "fasta", "fastq"),
  append = FALSE,
  compress = endsWith(tolower(outfile), ".gz")
) {
  resolved <- resolve_extract_index(index, file)
  index <- resolved$index
  file <- resolved$file
  index <- validate_index(index)
  live_index <- ensure_live_index_ptr(index)
  index <- live_index$index
  index_ptr <- live_index$ptr

  if (!is.character(outfile) || length(outfile) != 1L || !nzchar(outfile)) {
    stop("`outfile` must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.logical(append) || length(append) != 1L || is.na(append)) {
    stop("`append` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(compress) || length(compress) != 1L || is.na(compress)) {
    stop("`compress` must be TRUE or FALSE.", call. = FALSE)
  }

  if (length(seq_idx) < 1L) {
    return(invisible(normalizePath(outfile, winslash = "/", mustWork = FALSE)))
  }
  if (!is.numeric(seq_idx)) {
    stop("`seq_idx` must be numeric/integer-like.", call. = FALSE)
  }
  if (any(is.na(seq_idx) | seq_idx < 1 | seq_idx != floor(seq_idx))) {
    stop(
      "`seq_idx` must contain positive whole numbers (1-based).",
      call. = FALSE
    )
  }

  n_records <- as.numeric(index$n_records)
  if (any(seq_idx > n_records)) {
    stop(
      sprintf("Some seq_idx exceed available records (%s).", n_records),
      call. = FALSE
    )
  }

  files <- if (is.null(file)) {
    as.character(index$files)
  } else {
    files <- validate_input_files(file)
    if (length(files) != length(index$files)) {
      stop(
        paste0(
          "Override `file` must contain the same number of files as ",
          "index$files."
        ),
        call. = FALSE
      )
    }
    files
  }

  type <- match.arg(type)
  # nolint nextline: object_usage_linter
  cpp_extract_sequences_to_file(
    files = files,
    source_type = index$format,
    ids_zero_based = as.numeric(seq_idx - 1),
    index_ptr_sexp = index_ptr,
    output_type = type,
    outfile = outfile,
    append = append,
    compress = compress
  )

  invisible(normalizePath(outfile, winslash = "/", mustWork = FALSE))
}
