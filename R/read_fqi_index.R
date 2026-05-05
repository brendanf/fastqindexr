#' @noRd
normalize_existing_or_keep <- function(paths) {
  paths <- as.character(paths)
  existing <- file.exists(paths)
  if (any(existing)) {
    paths[existing] <- normalizePath(paths[existing], mustWork = TRUE)
  }
  paths
}

#' @noRd
deduce_files_from_fqi <- function(fqi_path) {
  fqi_path <- as.character(fqi_path)
  has_suffix <- grepl("\\.fqi$", fqi_path, ignore.case = TRUE)
  deduced <- ifelse(
    has_suffix,
    sub("\\.fqi$", "", fqi_path, ignore.case = TRUE),
    paste0(fqi_path, ".gz")
  )
  deduced <- normalize_existing_or_keep(deduced)
  missing <- deduced[!file.exists(deduced)]
  if (length(missing) > 0L) {
    warning(
      paste0(
        "Deduced source file(s) not found: ",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  deduced
}

#' @noRd
detect_type_from_gz <- function(path) {
  con <- try(gzfile(path, "rt"), silent = TRUE)
  if (inherits(con, "try-error")) {
    return(NA_character_)
  }
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) {
      return(NA_character_)
    }
    if (nzchar(line)) {
      first <- substr(line, 1L, 1L)
      if (identical(first, "@")) {
        return("fastq")
      }
      if (identical(first, ">")) {
        return("fasta")
      }
      return(NA_character_)
    }
  }
}

#' @noRd
resolve_type_for_files <- function(type, files) {
  if (!identical(type, "auto")) {
    return(type)
  }

  existing <- file.exists(files)
  detected <- rep(NA_character_, length(files))
  if (any(existing)) {
    detected[existing] <- vapply(
      files[existing],
      detect_type_from_gz,
      FUN.VALUE = character(1),
      USE.NAMES = FALSE
    )
  }

  unknown <- files[is.na(detected)]
  if (length(unknown) > 0L) {
    warning(
      paste0(
        "Could not determine type from file(s): ",
        paste(unknown, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  detected_unique <- unique(stats::na.omit(detected))
  if (length(detected_unique) > 1L) {
    stop(
      "Type auto-detection found mixed FASTA/FASTQ inputs across files.",
      call. = FALSE
    )
  }

  if (length(detected_unique) == 0L) {
    warning(
      "Type could not be determined; defaulting to 'fastq'.",
      call. = FALSE
    )
    return("fastq")
  }

  detected_unique[[1L]]
}

#' Read one or more FastqIndEx `.fqi` binary indexes
#'
#' Reads one or more `.fqi` index files produced by the FastqIndEx CLI and
#' converts them into a `fastqindexr_index` object that works with
#' [extract_sequences()].
#'
#' @param fqi_path Character vector of `.fqi` file paths.
#' @param files Optional character vector of indexed gz file paths. When omitted,
#'   paths are deduced using the convention that `.fqi` is appended to the source
#'   filename.
#' @param type Record format: `"auto"`, `"fasta"`, or `"fastq"`. For `"auto"`,
#'   type is inferred from the first non-empty line of each available source
#'   file (`>` for FASTA, `@` for FASTQ). If detection fails for all files, a
#'   warning is issued and `"fastq"` is used.
#'
#' @return A `fastqindexr_index` object with subclass
#'   `fastqindexr_gzip_index`.
#'
#' @export
read_fqi_index <- function(
  fqi_path,
  files = NULL,
  type = c("auto", "fasta", "fastq")
) {
  type <- match.arg(type)
  fqi_files <- validate_input_files(fqi_path)

  if (is.null(files)) {
    files <- deduce_files_from_fqi(fqi_files)
  } else {
    files <- validate_input_files(files)
  }

  if (length(files) != length(fqi_files)) {
    stop("`files` must have the same length as `fqi_path`.", call. = FALSE)
  }

  resolved_type <- resolve_type_for_files(type, files)
  index <- cpp_read_fqi_index(
    fqi_files = fqi_files,
    files = as.character(files),
    type = resolved_type
  )
  cache <- new.env(parent = emptyenv())
  cache$index_ptr <- index$index_ptr
  index$index_ptr <- NULL
  index$._cache <- cache
  assign_index_classes(index, default_subclass = "fastqindexr_gzip_index")
}
