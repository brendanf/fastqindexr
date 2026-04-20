#' @noRd
validate_input_files <- function(files) {
  if (length(files) < 1L) {
    stop("`files` must contain at least one path.", call. = FALSE)
  }
  files <- as.character(files)
  missing_files <- !file.exists(files)
  if (any(missing_files)) {
    missing_paths <- paste(files[missing_files], collapse = ", ")
    stop(
      sprintf("Missing file(s): %s", missing_paths),
      call. = FALSE
    )
  }
  normalizePath(files, mustWork = TRUE)
}

#' Build an in-memory index for gzipped FASTA or FASTQ
#'
#' Scans one or more gzipped files and builds a random-access style index. When
#' multiple files are given, they are treated as one logical concatenated
#' stream: record IDs follow file order (first file first, then the next, and
#' so on).
#'
#' @param files Character vector of paths to existing
#'   **gzip-compressed** (`.gz`)
#'   FASTA or FASTQ files. Paths are normalized to absolute paths.
#' @param type Format of the records: `"auto"` (infer from the first non-empty
#'   line of the first file: `>` for FASTA, `@` for FASTQ), `"fasta"`, or
#'   `"fastq"`.
#' @param index_stride_bytes Reserved for future use; currently ignored.
#'
#' @return An object of class `fastqindexr_index` (a list with additional
#'   attributes). Useful components include:
#' \describe{
#'   \item{format}{`"fasta"` or `"fastq"`.}
#'   \item{files}{Character vector of indexed paths.}
#'   \item{n_records}{Total number of records across all files.}
#'   \item{file_record_offsets}{Cumulative record offsets (length
#'     `length(files) + 1`); global record `i` maps to a file using these
#'     boundaries.}
#'   \item{file_record_counts}{Records per file.}
#'   \item{file_info}{Data frame with `path`, `size`, and `mtime` at index
#'     time.}
#'   \item{index_token}{Internal handle used with [extract_sequences()].}
#'   \item{record_size}{Lines per record (2 for FASTA, 4 for FASTQ).}
#' }
#'
#' @section FASTA limitation:
#' Extraction assumes **one sequence line per record** (header + sequence line).
#'
#' @seealso [extract_sequences()] to retrieve sequences by 1-based record ID.
#'
#' @examples
#' path <- tempfile(fileext = ".fasta.gz")
#' con <- gzfile(path, "wt")
#' writeLines(c(">s1", "AA", ">s2", "CC", ">s3", "GG"), con)
#' close(con)
#' idx <- create_index(path, type = "fasta")
#' extract_sequences(idx, seq_idx = c(3, 1))
#' unlink(path)
#'
#' @export
create_index <- function(
  files,
  type = c("auto", "fasta", "fastq"),
  index_stride_bytes = NULL
) {
  type <- match.arg(type)
  files <- validate_input_files(files)

  if (!is.null(index_stride_bytes)) {
    is_valid_stride <- is.numeric(index_stride_bytes) &&
      length(index_stride_bytes) == 1L &&
      !is.na(index_stride_bytes)
    if (!is_valid_stride) {
      stop(
        "`index_stride_bytes` must be NULL or a single numeric value.",
        call. = FALSE
      )
    }
  }

  index <- cpp_create_index(files = files, type = type)
  class(index) <- c("fastqindexr_index", class(index))
  index
}

#' @rdname create_index
#' @param x A `fastqindexr_index` object from [create_index()].
#' @param ... Ignored.
#' @export
print.fastqindexr_index <- function(x, ...) {
  cat("<fastqindexr_index>\n")
  cat("  format:", x$format, "\n")
  cat("  files:", length(x$files), "\n")
  cat("  records:", format(x$n_records, scientific = FALSE), "\n")
  invisible(x)
}
