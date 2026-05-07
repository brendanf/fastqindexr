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

#' @noRd
assign_index_classes <- function(index, default_subclass = NULL) {
  if (is.null(index$file_compression)) {
    compression <- character()
  } else {
    compression <- as.character(index$file_compression)
  }
  compression <- unique(compression)
  compression <- compression[!is.na(compression) & nzchar(compression)]

  subclass <- default_subclass
  if (length(compression) > 0L) {
    if (length(compression) != 1L) {
      stop(
        "Index payload contains mixed file compression values; use homogeneous inputs.",
        call. = FALSE
      )
    }
    subclass <- switch(
      compression[[1L]],
      gzip = "fastqindexr_gzip_index",
      plain = "fastqindexr_plain_index",
      NULL
    )
  }
  if (is.null(subclass)) {
    stop(
      "Could not determine index subclass from `file_compression`.",
      call. = FALSE
    )
  }

  class(index) <- c(subclass, "fastqindexr_index", "list")
  index
}

#' Build an in-memory index for FASTA or FASTQ
#'
#' Scans one or more files and builds a random-access style index. Gzip inputs
#' use block-style indexing; uncompressed inputs use byte offsets for record
#' starts. When multiple files are given, they are treated as one logical
#' concatenated stream: record IDs follow file order (first file first, then
#' the next, and so on).
#'
#' @param files Character vector of paths to existing FASTA or FASTQ files
#'   (gzip-compressed or plain). Paths are normalized to absolute paths. A single
#'   call must use one compression type only (all gzip or all plain).
#' @param type Format of the records: `"auto"` (infer from the first non-empty
#'   line of the first file: `>` for FASTA, `@` for FASTQ), `"fasta"`, or
#'   `"fastq"`.
#' @param index_stride_bytes Reserved for future use; currently ignored.
#'
#' @return A `fastqindexr_index` object (list-based), with subclass
#'   `fastqindexr_gzip_index` or `fastqindexr_plain_index`. Useful components
#'   include:
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
#'   \item{index_payload}{Serialized internal index payload used to restore
#'     native state after serialization/deserialization.}
#'   \item{record_size}{Internal extraction stride (2 for FASTA, 4 for FASTQ;
#'     FASTA logical records are header-driven and may span multiple lines).}
#'   \item{file_compression}{Per-file `"gzip"` or `"plain"`.}
#' }
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
  cache <- new.env(parent = emptyenv())
  cache$index_ptr <- index$index_ptr
  index$index_ptr <- NULL
  index$._cache <- cache
  assign_index_classes(index)
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
