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
#' Returns rows in **the same order as `ids`**, including duplicate IDs. IDs are
#' **1-based** record indices in the logical concatenated stream defined when
#' the index was built (see [create_index()]).
#'
#' @param index A `fastqindexr_index` object from [create_index()].
#' @param seq_idx Numeric vector of record IDs (positive whole numbers). Values
#'   are coerced via [as.numeric()]; `NA` and values outside `1` … `n_records`
#    are errors.
#' @param file Optional character vector of file paths overriding those
#'   stored in
#'   `index$files`. Must be the same length as `index$files` when provided; each
#'   path must exist. Use this when data moved after indexing.
#'
#' @return A data frame:
#' \describe{
#'   \item{FASTA}{Columns `seq_id` and `seq`.}
#'   \item{FASTQ}{Columns `seq_id`, `seq`, and `qual`.}
#' }
#' For `length(ids) == 0`, an empty data frame with the appropriate columns is
#' returned.
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
extract_sequences <- function(index, seq_idx, file = NULL) {
  index <- validate_index(index)
  live_index <- ensure_live_index_ptr(index)
  index <- live_index$index
  index_ptr <- live_index$ptr

  if (length(seq_idx) < 1L) {
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

  # nolint nextline: object_usage_linter
  result <- cpp_extract_sequences(
    files = files,
    type = index$format,
    ids_zero_based = as.numeric(seq_idx - 1),
    index_ptr_sexp = index_ptr
  )

  if (identical(index$format, "fastq")) {
    out <- data.frame(
      seq_id = as.character(result$seq_id),
      seq = as.character(result$seq),
      qual = as.character(result$qual),
      stringsAsFactors = FALSE
    )
  } else {
    out <- data.frame(
      seq_id = as.character(result$seq_id),
      seq = as.character(result$seq),
      stringsAsFactors = FALSE
    )
  }
  out
}
