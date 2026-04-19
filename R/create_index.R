validate_input_files <- function(files) {
  if (length(files) < 1L) {
    stop("`files` must contain at least one path.", call. = FALSE)
  }
  files <- as.character(files)
  missing_files <- !file.exists(files)
  if (any(missing_files)) {
    stop(
      sprintf("Missing file(s): %s", paste(files[missing_files], collapse = ", ")),
      call. = FALSE
    )
  }
  normalizePath(files, mustWork = TRUE)
}

create_index <- function(files, type = c("auto", "fasta", "fastq"), index_stride_bytes = NULL) {
  type <- match.arg(type)
  files <- validate_input_files(files)

  if (!is.null(index_stride_bytes)) {
    if (!is.numeric(index_stride_bytes) || length(index_stride_bytes) != 1L || is.na(index_stride_bytes)) {
      stop("`index_stride_bytes` must be NULL or a single numeric value.", call. = FALSE)
    }
  }

  index <- cpp_create_index(files = files, type = type)
  class(index) <- c("fastqindexr_index", class(index))
  index
}

print.fastqindexr_index <- function(x, ...) {
  cat("<fastqindexr_index>\n")
  cat("  format:", x$format, "\n")
  cat("  files:", length(x$files), "\n")
  cat("  records:", format(x$n_records, scientific = FALSE), "\n")
  invisible(x)
}
