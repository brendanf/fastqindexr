#' @noRd
validate_index <- function(index) {
  if (!inherits(index, "fastqindexr_index")) {
    stop("`index` must be a fastqindexr_index object.", call. = FALSE)
  }
  required <- c("files", "format", "n_records", "file_record_offsets", "index_token")
  missing <- setdiff(required, names(index))
  if (length(missing) > 0L) {
    stop(
      sprintf("Malformed index object. Missing: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  index
}

#' Extract sequence records by ID from indexed gzipped FASTA or FASTQ
#'
#' Returns rows in **the same order as `ids`**, including duplicate IDs. IDs are
#' **1-based** record indices in the logical concatenated stream defined when
#' the index was built (see [create_index()]).
#'
#' @param index A `fastqindexr_index` object from [create_index()].
#' @param ids Numeric vector of record IDs (positive whole numbers). Values are
#'   coerced via [as.numeric()]; `NA` and values outside `1` … `n_records` are
#'   errors.
#' @param file Optional character vector of file paths overriding those stored in
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
#' extract_sequences(idx, ids = c(2, 2, 1))
#' unlink(path)
#'
#' @export
extract_sequences <- function(index, ids, file = NULL) {
  index <- validate_index(index)

  if (length(ids) < 1L) {
    if (identical(index$format, "fastq")) {
      return(data.frame(seq_id = character(), seq = character(), qual = character(), stringsAsFactors = FALSE))
    }
    return(data.frame(seq_id = character(), seq = character(), stringsAsFactors = FALSE))
  }

  if (!is.numeric(ids)) {
    stop("`ids` must be numeric/integer-like.", call. = FALSE)
  }
  if (any(is.na(ids) | ids < 1 | ids != floor(ids))) {
    stop("`ids` must contain positive whole numbers (1-based).", call. = FALSE)
  }

  n_records <- as.numeric(index$n_records)
  if (any(ids > n_records)) {
    stop(sprintf("Some ids exceed available records (%s).", n_records), call. = FALSE)
  }

  files <- if (is.null(file)) {
    as.character(index$files)
  } else {
    files <- validate_input_files(file)
    if (length(files) != length(index$files)) {
      stop("Override `file` must contain the same number of files as index$files.", call. = FALSE)
    }
    files
  }

  result <- cpp_extract_sequences(
    files = files,
    type = index$format,
    ids_zero_based = as.numeric(ids - 1),
    index_token = index$index_token
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
