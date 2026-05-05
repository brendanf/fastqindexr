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
  compression <- as.character(index$file_compression)
  compression <- unique(compression[!is.na(compression) & nzchar(compression)])
  if (
    length(compression) != 1L || !(compression[[1L]] %in% c("gzip", "plain"))
  ) {
    stop(
      paste0(
        "Malformed index object. `file_compression` must be one homogeneous ",
        "value in {\"gzip\", \"plain\"}."
      ),
      call. = FALSE
    )
  }
  if (
    inherits(index, "fastqindexr_gzip_index") &&
      !identical(compression[[1L]], "gzip")
  ) {
    stop(
      "Malformed index object. gzip subclass requires gzip compression.",
      call. = FALSE
    )
  }
  if (
    inherits(index, "fastqindexr_plain_index") &&
      !identical(compression[[1L]], "plain")
  ) {
    stop(
      "Malformed index object. plain subclass requires plain compression.",
      call. = FALSE
    )
  }
  index
}

#' @noRd
resolve_extract_index <- function(index, file) {
  if (is.null(index)) {
    return(list(index = NULL, file = file))
  }
  if (inherits(index, "fastqindexr_index")) {
    return(list(index = index, file = file))
  }
  if (is.character(index)) {
    loaded <- read_fqi_index(fqi_path = index, files = file, type = "auto")
    return(list(index = loaded, file = NULL))
  }
  stop(
    "`index` must be NULL, a fastqindexr_index object, or path(s) to `.fqi` file(s).",
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

#' @noRd
validate_seq_idx <- function(seq_idx, n_records, arg_name = "`seq_idx`") {
  if (!is.numeric(seq_idx)) {
    stop(sprintf("%s must be numeric/integer-like.", arg_name), call. = FALSE)
  }
  if (any(is.na(seq_idx) | seq_idx < 1 | seq_idx != floor(seq_idx))) {
    stop(
      sprintf("%s must contain positive whole numbers (1-based).", arg_name),
      call. = FALSE
    )
  }
  if (any(seq_idx > n_records)) {
    stop(
      sprintf("Some %s exceed available records (%s).", arg_name, n_records),
      call. = FALSE
    )
  }
  as.numeric(seq_idx)
}

#' @noRd
create_empty_extract_file <- function(path, compress) {
  if (compress) {
    con <- gzfile(path, open = "wb")
  } else {
    con <- file(path, open = "wb")
  }
  close(con)
}

#' @noRd
write_sequences_df_to_file <- function(
  df,
  outfile,
  out_type,
  compress,
  append
) {
  mode <- if (append) "ab" else "wb"
  con <- if (compress) {
    gzfile(outfile, open = mode)
  } else {
    file(outfile, open = mode)
  }
  on.exit(close(con), add = TRUE)
  if (identical(out_type, "fastq")) {
    for (i in seq_len(nrow(df))) {
      cat(
        paste0(
          "@",
          df$seq_id[[i]],
          "\n",
          df$seq[[i]],
          "\n+\n",
          df$qual[[i]],
          "\n"
        ),
        file = con
      )
    }
  } else {
    for (i in seq_len(nrow(df))) {
      cat(
        paste0(">", df$seq_id[[i]], "\n", df$seq[[i]], "\n"),
        file = con
      )
    }
  }
  invisible(outfile)
}

#' @noRd
resolve_region_merge_tuning <- function() {
  bridge_gap <- getOption("fastqindexr.max_bridge_gap", 512)
  region_records <- getOption("fastqindexr.max_region_records", 100000)
  extract_mode <- getOption("fastqindexr.extract_mode", "indexed")
  diagnostics <- getOption("fastqindexr.extract_diagnostics", FALSE)

  if (
    !is.numeric(bridge_gap) ||
      length(bridge_gap) != 1L ||
      is.na(bridge_gap) ||
      !is.finite(bridge_gap) ||
      bridge_gap < 0 ||
      bridge_gap != floor(bridge_gap)
  ) {
    stop(
      paste0(
        "Option `fastqindexr.max_bridge_gap` must be a finite ",
        "non-negative whole number."
      ),
      call. = FALSE
    )
  }
  if (
    !is.numeric(region_records) ||
      length(region_records) != 1L ||
      is.na(region_records) ||
      !is.finite(region_records) ||
      region_records <= 0 ||
      region_records != floor(region_records)
  ) {
    stop(
      paste0(
        "Option `fastqindexr.max_region_records` must be a finite ",
        "positive whole number."
      ),
      call. = FALSE
    )
  }
  if (
    !is.character(extract_mode) ||
      length(extract_mode) != 1L ||
      is.na(extract_mode) ||
      !(extract_mode %in% c("indexed", "sequential_only", "sequential"))
  ) {
    stop(
      paste0(
        "Option `fastqindexr.extract_mode` must be one of ",
        "`\"indexed\"`, `\"sequential\"`, or `\"sequential_only\"` (deprecated)."
      ),
      call. = FALSE
    )
  }
  if (identical(extract_mode, "sequential_only")) {
    warning(
      "`options(fastqindexr.extract_mode = \"sequential_only\")` is deprecated; use `\"sequential\"`.",
      call. = FALSE
    )
    extract_mode <- "sequential"
  }
  if (
    !is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)
  ) {
    stop(
      "Option `fastqindexr.extract_diagnostics` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  list(
    max_bridge_gap = as.numeric(bridge_gap),
    max_region_records = as.numeric(region_records),
    extract_mode = extract_mode,
    diagnostics = diagnostics
  )
}

#' @noRd
resolve_source_type_no_index <- function(files, type) {
  type <- match.arg(type, c("auto", "fasta", "fastq"))
  if (identical(type, "auto")) {
    cpp_detect_input_format(files[[1L]])
  } else {
    type
  }
}

#' @noRd
streaming_record_offsets <- function(index, files, source_type) {
  if (!is.null(index)) {
    return(as.numeric(index$file_record_offsets))
  }
  cpp_scan_record_offsets(files = files, type = source_type)
}

#' @noRd
extraction_uses_streaming <- function(mode, has_index, option_extract_mode) {
  if (identical(mode, "sequential")) {
    return(TRUE)
  }
  if (identical(mode, "indexed")) {
    return(FALSE)
  }
  if (!has_index) {
    return(TRUE)
  }
  identical(option_extract_mode, "sequential")
}

#' @noRd
apply_extract_renumber <- function(seq_ids, renumber) {
  renumber <- match.arg(renumber, c("none", "zero_based", "one_based"))
  if (identical(renumber, "none")) {
    return(seq_ids)
  }
  n <- seq_along(seq_ids)
  if (identical(renumber, "zero_based")) {
    return(as.character(n - 1L))
  }
  as.character(n)
}

#' Partition sequence IDs into stable batches
#'
#' Splits a numeric vector of 1-based sequence IDs into `n_parts` partitions
#' while preserving order within each partition. Duplicate IDs are preserved.
#'
#' @param seq_idx Numeric vector of positive whole-number sequence IDs.
#' @param n_parts Number of partitions to create (positive whole number).
#' @param strategy Partitioning strategy: `"contiguous"` or `"round_robin"`.
#' @param drop_empty If `TRUE`, empty partitions are removed from the output.
#'
#' @return A list of integer vectors containing partitioned IDs.
#'
#' @export
partition_seq_idx <- function(
  seq_idx,
  n_parts,
  strategy = c("contiguous", "round_robin"),
  drop_empty = TRUE
) {
  strategy <- match.arg(strategy)
  if (!is.numeric(seq_idx)) {
    stop("`seq_idx` must be numeric/integer-like.", call. = FALSE)
  }
  if (any(is.na(seq_idx) | seq_idx < 1 | seq_idx != floor(seq_idx))) {
    stop(
      "`seq_idx` must contain positive whole numbers (1-based).",
      call. = FALSE
    )
  }
  if (
    !is.numeric(n_parts) ||
      length(n_parts) != 1L ||
      is.na(n_parts) ||
      n_parts < 1 ||
      n_parts != floor(n_parts)
  ) {
    stop("`n_parts` must be a positive whole number.", call. = FALSE)
  }
  if (
    !is.logical(drop_empty) || length(drop_empty) != 1L || is.na(drop_empty)
  ) {
    stop("`drop_empty` must be TRUE or FALSE.", call. = FALSE)
  }

  seq_idx <- as.integer(seq_idx)
  n_parts <- as.integer(n_parts)
  out <- replicate(n_parts, integer(), simplify = FALSE)
  if (length(seq_idx) > 0L) {
    if (identical(strategy, "contiguous")) {
      base_size <- length(seq_idx) %/% n_parts
      n_larger <- length(seq_idx) %% n_parts
      start <- 1L
      for (i in seq_len(n_parts)) {
        size_i <- base_size + if (i <= n_larger) 1L else 0L
        if (size_i > 0L) {
          end <- start + size_i - 1L
          out[[i]] <- seq_idx[start:end]
          start <- end + 1L
        }
      }
    } else {
      slot <- ((seq_along(seq_idx) - 1L) %% n_parts) + 1L
      out <- split(seq_idx, slot)
      out <- out[as.character(seq_len(n_parts))]
    }
  }

  if (drop_empty) {
    out[vapply(out, length, integer(1L)) > 0L]
  } else {
    out
  }
}

#' Extract sequence records by ID from indexed FASTA or FASTQ
#'
#' Returns rows in **the same order as `seq_idx`**, including duplicate IDs. IDs
#' are
#' **1-based** record indices in the logical concatenated stream defined when
#' the index was built (see [create_index()]).
#'
#' @param index Either a `fastqindexr_index` object (from [create_index()] or
#'   [read_fqi_index()]) or one or more `.fqi` file paths.
#' @param seq_idx Numeric vector of record IDs (positive whole numbers), or a
#'   **list** of such vectors (one result per list element, same order). Values
#'   are coerced via [as.numeric()]; `NA` and values outside `1` &hellip; `n_records`
#'   are errors.
#' @param file When `index` is non-`NULL`, optional path override(s) for
#'   `index$files`. When `index` is `NULL`, **required** source file path(s) for
#'   streaming extraction (`mode` `"auto"` or `"sequential"`).
#'   For `.fqi` path input, this provides indexed file path(s) passed to
#'   [read_fqi_index()]. If omitted for `.fqi` input, file paths are deduced from
#'   the `.fqi` names.
#' @param return In-memory return shape: `"data.frame"` (the default), `"list"`, or
#'   `"seq"`. With `"data.frame"`, the result is as before. With `"list"`, returns
#'   a list with `seq_id` and `seq` (and `qual` for FASTQ). With `"seq"`, returns
#'   a character vector of sequences with `seq_id` as names, preserving order
#'   and allowing duplicate names. For indexed FASTQ, `return = "seq"` skips
#'   reading quality lines from the source.
#' @param mode `"auto"` uses index-backed extraction when `index` is non-`NULL`
#'   and `getOption("fastqindexr.extract_mode", "indexed")` is `"indexed"`;
#'   otherwise (no index, or option `"sequential"`) reads by streaming scans.
#'   `"indexed"` requires `index` and offset-based I/O. `"sequential"` always
#'   streams and does not require a live native index pointer (works with
#'   `index` fields after [readRDS()] if paths remain valid).
#' @param type When `index` is `NULL`, passed to format detection: `"auto"`,
#'   `"fasta"`, or `"fastq"`.
#' @param renumber Renames emitted `seq_id` / FASTA headers: `"none"`,
#'   `"zero_based"`, or `"one_based"` (output order positions).
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
#' If `seq_idx` is a list, the return value is a list of results (same length
#' and order as `seq_idx`), each obeying `return` and `renumber` for that
#' partition.
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
  index = NULL,
  seq_idx,
  file = NULL,
  return = c("data.frame", "list", "seq"),
  mode = c("auto", "indexed", "sequential"),
  type = c("auto", "fasta", "fastq"),
  renumber = c("none", "zero_based", "one_based")
) {
  return <- match.arg(return)
  mode <- match.arg(mode)
  renumber <- match.arg(renumber)
  tuning <- resolve_region_merge_tuning()
  resolved <- resolve_extract_index(index, file)
  index <- resolved$index
  file <- resolved$file
  has_index <- !is.null(index)
  if (identical(mode, "indexed") && !has_index) {
    stop("`mode = \"indexed\"` requires a non-NULL `index`.", call. = FALSE)
  }
  if (has_index) {
    index <- validate_index(index)
  }

  empty_out_df <- function(fmt) {
    if (identical(fmt, "fastq")) {
      return(data.frame(
        seq_id = character(),
        seq = character(),
        qual = character(),
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      seq_id = character(),
      seq = character(),
      stringsAsFactors = FALSE
    )
  }
  empty_out_list <- function(fmt) {
    if (identical(fmt, "fastq")) {
      return(list(
        seq_id = character(),
        seq = character(),
        qual = character()
      ))
    }
    list(seq_id = character(), seq = character())
  }

  if (is.list(seq_idx)) {
    out <- vector("list", length(seq_idx))
    for (i in seq_along(seq_idx)) {
      out[[i]] <- extract_sequences(
        index = index,
        seq_idx = seq_idx[[i]],
        file = file,
        return = return,
        mode = mode,
        type = type,
        renumber = renumber
      )
    }
    return(out)
  }

  source_type <- if (has_index) {
    index$format
  } else {
    files_pre <- validate_input_files(file)
    resolve_source_type_no_index(files_pre, type)
  }

  use_streaming <- extraction_uses_streaming(
    mode,
    has_index,
    tuning$extract_mode
  )

  if (length(seq_idx) < 1L) {
    if (return == "list") {
      return(empty_out_list(source_type))
    }
    if (return == "seq") {
      return(structure(character(0L), names = character(0L)))
    }
    return(empty_out_df(source_type))
  }

  files <- if (has_index) {
    if (is.null(file)) {
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
  } else {
    validate_input_files(file)
  }

  n_records <- if (has_index) {
    as.numeric(index$n_records)
  } else {
    roff <- streaming_record_offsets(NULL, files, source_type)
    as.numeric(roff[[length(roff)]])
  }
  seq_idx <- validate_seq_idx(seq_idx = seq_idx, n_records = n_records)

  is_fastq <- identical(source_type, "fastq")
  include_qual <- is_fastq && (return %in% c("data.frame", "list"))

  if (use_streaming) {
    offsets <- streaming_record_offsets(index, files, source_type)
    # nolint nextline: object_usage_linter
    result <- cpp_extract_sequences_streaming(
      files = files,
      type = source_type,
      ids_zero_based = as.numeric(seq_idx - 1),
      record_offsets_r = offsets,
      include_qual = include_qual,
      diagnostics = tuning$diagnostics
    )
  } else {
    live_index <- ensure_live_index_ptr(index)
    index <- live_index$index
    index_ptr <- live_index$ptr
    # nolint nextline: object_usage_linter
    result <- cpp_extract_sequences(
      files = files,
      type = index$format,
      ids_zero_based = as.numeric(seq_idx - 1),
      index_ptr_sexp = index_ptr,
      include_qual = include_qual,
      max_bridge_gap = tuning$max_bridge_gap,
      max_region_records = tuning$max_region_records,
      extract_mode = tuning$extract_mode,
      diagnostics = tuning$diagnostics
    )
  }
  diag <- attr(result, "fastqindexr_diagnostics", exact = TRUE)

  if (!identical(renumber, "none")) {
    new_ids <- apply_extract_renumber(
      as.character(result$seq_id),
      renumber
    )
    result$seq_id <- new_ids
  }

  if (return == "data.frame") {
    if (is_fastq) {
      return(
        structure(
          data.frame(
            seq_id = result$seq_id,
            seq = result$seq,
            qual = result$qual,
            stringsAsFactors = FALSE
          ),
          fastqindexr_diagnostics = diag
        )
      )
    }
    return(
      structure(
        data.frame(
          seq_id = result$seq_id,
          seq = result$seq,
          stringsAsFactors = FALSE
        ),
        fastqindexr_diagnostics = diag
      )
    )
  }
  if (return == "list") {
    if (is_fastq) {
      return(
        structure(
          list(
            seq_id = result$seq_id,
            seq = result$seq,
            qual = result$qual
          ),
          fastqindexr_diagnostics = diag
        )
      )
    }
    return(
      structure(
        list(seq_id = result$seq_id, seq = result$seq),
        fastqindexr_diagnostics = diag
      )
    )
  }
  structure(
    stats::setNames(result$seq, result$seq_id),
    fastqindexr_diagnostics = diag
  )
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
  index = NULL,
  seq_idx,
  file = NULL,
  outfile,
  type = c("auto", "fasta", "fastq"),
  append = FALSE,
  compress = endsWith(tolower(outfile), ".gz"),
  mode = c("auto", "indexed", "sequential"),
  input_type = c("auto", "fasta", "fastq"),
  renumber = c("none", "zero_based", "one_based")
) {
  mode <- match.arg(mode)
  renumber <- match.arg(renumber)
  tuning <- resolve_region_merge_tuning()
  resolved <- resolve_extract_index(index, file)
  index <- resolved$index
  file <- resolved$file
  has_index <- !is.null(index)
  if (identical(mode, "indexed") && !has_index) {
    stop("`mode = \"indexed\"` requires a non-NULL `index`.", call. = FALSE)
  }
  if (has_index) {
    index <- validate_index(index)
  }

  if (!is.logical(append) || length(append) != 1L || is.na(append)) {
    stop("`append` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(outfile) || length(outfile) < 1L || any(!nzchar(outfile))) {
    stop("`outfile` must be a non-empty character vector.", call. = FALSE)
  }
  is_partition_mode <- is.list(seq_idx) || length(outfile) > 1L

  source_type <- if (has_index) {
    index$format
  } else {
    files_pre <- validate_input_files(file)
    resolve_source_type_no_index(files_pre, input_type)
  }

  use_streaming <- extraction_uses_streaming(
    mode,
    has_index,
    tuning$extract_mode
  )

  files <- if (has_index) {
    if (is.null(file)) {
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
  } else {
    validate_input_files(file)
  }

  n_records <- if (has_index) {
    as.numeric(index$n_records)
  } else {
    roff <- streaming_record_offsets(NULL, files, source_type)
    as.numeric(roff[[length(roff)]])
  }

  type <- match.arg(type)

  if (!identical(renumber, "none")) {
    if (is_partition_mode) {
      if (!is.list(seq_idx)) {
        stop("Partitioned mode requires `seq_idx` to be a list.", call. = FALSE)
      }
      if (length(seq_idx) != length(outfile)) {
        stop("`outfile` must match `seq_idx` list length.", call. = FALSE)
      }
      compress_vec <- if (length(compress) == 1L) {
        rep(compress, length(outfile))
      } else {
        compress
      }
      out_norm <- normalizePath(outfile, winslash = "/", mustWork = FALSE)
      for (i in seq_along(seq_idx)) {
        ids_i <- seq_idx[[i]]
        if (length(ids_i) < 1L) {
          if (!append) {
            create_empty_extract_file(outfile[i], compress = compress_vec[i])
          }
          next
        }
        validate_seq_idx(ids_i, n_records, "`seq_idx` partition")
        df <- extract_sequences(
          index = index,
          seq_idx = ids_i,
          file = file,
          return = "data.frame",
          mode = mode,
          type = input_type,
          renumber = renumber
        )
        resolved_type <- if (type == "auto") source_type else type
        if (source_type == "fasta" && resolved_type == "fastq") {
          stop("Cannot emit FASTQ output from FASTA input.", call. = FALSE)
        }
        write_sequences_df_to_file(
          df,
          outfile[i],
          resolved_type,
          compress_vec[i],
          append
        )
      }
      return(invisible(out_norm))
    }
    if (length(seq_idx) < 1L) {
      return(invisible(normalizePath(
        outfile,
        winslash = "/",
        mustWork = FALSE
      )))
    }
    validate_seq_idx(seq_idx, n_records)
    df <- extract_sequences(
      index = index,
      seq_idx = seq_idx,
      file = file,
      return = "data.frame",
      mode = mode,
      type = input_type,
      renumber = renumber
    )
    resolved_type <- if (type == "auto") source_type else type
    if (source_type == "fasta" && resolved_type == "fastq") {
      stop("Cannot emit FASTQ output from FASTA input.", call. = FALSE)
    }
    write_sequences_df_to_file(df, outfile, resolved_type, compress, append)
    return(invisible(normalizePath(outfile, winslash = "/", mustWork = FALSE)))
  }

  live_index <- if (!use_streaming && has_index) {
    ensure_live_index_ptr(index)
  } else {
    NULL
  }
  index_ptr <- if (!is.null(live_index)) {
    live_index$ptr
  } else {
    NULL
  }
  if (!is.null(live_index)) {
    index <- live_index$index
  }

  if (is_partition_mode) {
    if (!is.list(seq_idx)) {
      stop("Partitioned mode requires `seq_idx` to be a list.", call. = FALSE)
    }
    if (length(seq_idx) != length(outfile)) {
      stop("`outfile` must match `seq_idx` list length.", call. = FALSE)
    }
    if (
      !is.logical(compress) ||
        any(is.na(compress)) ||
        !(length(compress) %in% c(1L, length(outfile)))
    ) {
      stop(
        paste0(
          "`compress` must be TRUE/FALSE scalar or vector matching ",
          "`outfile` length."
        ),
        call. = FALSE
      )
    }
    compress_vec <- if (length(compress) == 1L) {
      rep(compress, length(outfile))
    } else {
      compress
    }
    out_norm <- normalizePath(outfile, winslash = "/", mustWork = FALSE)
    validated_parts <- vector("list", length(seq_idx))
    all_strict_increasing <- TRUE
    any_ids <- FALSE
    for (i in seq_along(seq_idx)) {
      ids_i <- seq_idx[[i]]
      if (length(ids_i) < 1L) {
        validated_parts[[i]] <- numeric()
        if (!append) {
          create_empty_extract_file(outfile[i], compress = compress_vec[i])
        }
        next
      }
      ids_i <- validate_seq_idx(
        seq_idx = ids_i,
        n_records = n_records,
        arg_name = "`seq_idx` partition"
      )
      validated_parts[[i]] <- ids_i
      any_ids <- TRUE
      if (length(ids_i) > 1L && any(diff(ids_i) <= 0)) {
        all_strict_increasing <- FALSE
      }
    }

    # fastqindexr change: partition fast path merges unique IDs across
    # partitions to build indexed extraction regions once, then scatters.
    all_ids <- unlist(validated_parts, use.names = FALSE)
    no_cross_partition_duplicates <- length(all_ids) == length(unique(all_ids))
    can_merge_partitions <- !append &&
      any_ids &&
      all_strict_increasing &&
      no_cross_partition_duplicates

    if (can_merge_partitions) {
      resolved_type <- if (type == "auto") source_type else type
      if (source_type == "fasta" && resolved_type == "fastq") {
        stop("Cannot emit FASTQ output from FASTA index input.", call. = FALSE)
      }
      include_qual <- identical(source_type, "fastq") &&
        identical(resolved_type, "fastq")
      union_ids <- sort(all_ids)
      if (use_streaming) {
        offsets <- streaming_record_offsets(index, files, source_type)
        # nolint nextline: object_usage_linter
        merged <- cpp_extract_sequences_streaming(
          files = files,
          type = source_type,
          ids_zero_based = as.numeric(union_ids - 1),
          record_offsets_r = offsets,
          include_qual = include_qual,
          diagnostics = tuning$diagnostics
        )
      } else {
        # nolint nextline: object_usage_linter
        merged <- cpp_extract_sequences(
          files = files,
          type = index$format,
          ids_zero_based = as.numeric(union_ids - 1),
          index_ptr_sexp = index_ptr,
          include_qual = include_qual,
          max_bridge_gap = tuning$max_bridge_gap,
          max_region_records = tuning$max_region_records,
          extract_mode = tuning$extract_mode,
          diagnostics = tuning$diagnostics
        )
      }
      render_one <- function(seq_id, seq, qual) {
        if (resolved_type == "fastq") {
          paste0("@", seq_id, "\n", seq, "\n+\n", qual, "\n")
        } else {
          paste0(">", seq_id, "\n", seq, "\n")
        }
      }
      rendered <- vapply(
        seq_along(union_ids),
        function(j) {
          render_one(
            merged$seq_id[[j]],
            merged$seq[[j]],
            if (isTRUE(include_qual)) merged$qual[[j]] else ""
          )
        },
        FUN.VALUE = character(1L),
        USE.NAMES = FALSE
      )
      names(rendered) <- as.character(union_ids)

      for (i in seq_along(validated_parts)) {
        ids_i <- validated_parts[[i]]
        if (length(ids_i) < 1L) {
          next
        }
        mode <- if (append) "ab" else "wb"
        con <- if (compress_vec[i]) {
          gzfile(outfile[i], open = mode)
        } else {
          file(outfile[i], open = mode)
        }
        cat(rendered[as.character(ids_i)], file = con, sep = "")
        close(con)
      }
      if (isTRUE(tuning$diagnostics)) {
        attr(out_norm, "fastqindexr_diagnostics") <- attr(
          merged,
          "fastqindexr_diagnostics",
          exact = TRUE
        )
      }
      return(invisible(out_norm))
    }

    for (i in seq_along(validated_parts)) {
      ids_i <- validated_parts[[i]]
      if (length(ids_i) < 1L) {
        next
      }
      if (use_streaming) {
        offsets <- streaming_record_offsets(index, files, source_type)
        # nolint nextline: object_usage_linter
        res <- cpp_extract_sequences_to_file_streaming(
          files = files,
          source_type = source_type,
          ids_zero_based = as.numeric(ids_i - 1),
          record_offsets_r = offsets,
          output_type = type,
          outfile = outfile[i],
          append = append,
          compress = compress_vec[i],
          diagnostics = tuning$diagnostics
        )
      } else {
        # nolint nextline: object_usage_linter
        res <- cpp_extract_sequences_to_file(
          files = files,
          source_type = index$format,
          ids_zero_based = as.numeric(ids_i - 1),
          index_ptr_sexp = index_ptr,
          output_type = type,
          outfile = outfile[i],
          append = append,
          compress = compress_vec[i],
          max_bridge_gap = tuning$max_bridge_gap,
          max_region_records = tuning$max_region_records,
          extract_mode = tuning$extract_mode,
          diagnostics = tuning$diagnostics
        )
      }
      if (isTRUE(tuning$diagnostics)) {
        attr(out_norm, "fastqindexr_diagnostics") <- attr(
          res,
          "fastqindexr_diagnostics",
          exact = TRUE
        )
      }
    }
    return(invisible(out_norm))
  }

  if (!is.logical(compress) || length(compress) != 1L || is.na(compress)) {
    stop("`compress` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(seq_idx) < 1L) {
    return(invisible(normalizePath(outfile, winslash = "/", mustWork = FALSE)))
  }
  seq_idx <- validate_seq_idx(seq_idx = seq_idx, n_records = n_records)
  if (use_streaming) {
    offsets <- streaming_record_offsets(index, files, source_type)
    # nolint nextline: object_usage_linter
    res <- cpp_extract_sequences_to_file_streaming(
      files = files,
      source_type = source_type,
      ids_zero_based = as.numeric(seq_idx - 1),
      record_offsets_r = offsets,
      output_type = type,
      outfile = outfile,
      append = append,
      compress = compress,
      diagnostics = tuning$diagnostics
    )
  } else {
    # nolint nextline: object_usage_linter
    res <- cpp_extract_sequences_to_file(
      files = files,
      source_type = index$format,
      ids_zero_based = as.numeric(seq_idx - 1),
      index_ptr_sexp = index_ptr,
      output_type = type,
      outfile = outfile,
      append = append,
      compress = compress,
      max_bridge_gap = tuning$max_bridge_gap,
      max_region_records = tuning$max_region_records,
      extract_mode = tuning$extract_mode,
      diagnostics = tuning$diagnostics
    )
  }
  out_norm <- normalizePath(outfile, winslash = "/", mustWork = FALSE)
  if (isTRUE(tuning$diagnostics)) {
    attr(out_norm, "fastqindexr_diagnostics") <- attr(
      res,
      "fastqindexr_diagnostics",
      exact = TRUE
    )
  }

  invisible(out_norm)
}

#' Extract sequence records as a DNAStringSet
#'
#' Returns sequences as a `Biostrings::DNAStringSet` while preserving request
#' order and duplicate IDs. Extraction is processed in chunks to keep memory
#' usage bounded on large requests.
#'
#' @param index Either a `fastqindexr_index` object (from [create_index()] or
#'   [read_fqi_index()]) or one or more `.fqi` file paths.
#' @param seq_idx Numeric vector of record IDs (positive whole numbers).
#' @param file Optional character vector of file paths overriding those stored
#'   in `index$files` for object input.
#' @param renumber Name handling for returned sequences: `"none"` keeps source
#'   IDs, `"zero_based"` uses `"0"`, `"1"`, ..., and `"one_based"` uses
#'   `"1"`, `"2"`, ....
#' @param chunk_chars Approximate number of sequence characters to process per
#'   chunk when building the return object. Larger values may reduce overhead
#'   but increase temporary memory usage. Must be a finite positive number (not
#'   `NA`, not `Inf`).
#'
#' @return A `Biostrings::DNAStringSet`.
#'
#' @export
extract_sequences_dnastringset <- function(
  index = NULL,
  seq_idx,
  file = NULL,
  renumber = c("none", "zero_based", "one_based"),
  chunk_chars = 1e7,
  mode = c("auto", "indexed", "sequential"),
  type = c("auto", "fasta", "fastq")
) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package `Biostrings` is required.", call. = FALSE)
  }
  renumber <- match.arg(renumber)
  mode <- match.arg(mode)
  tuning <- resolve_region_merge_tuning()
  if (
    !is.numeric(chunk_chars) ||
      length(chunk_chars) != 1L ||
      is.na(chunk_chars) ||
      !is.finite(chunk_chars) ||
      chunk_chars <= 0
  ) {
    stop(
      "`chunk_chars` must be a finite positive numeric scalar.",
      call. = FALSE
    )
  }
  resolved <- resolve_extract_index(index, file)
  index <- resolved$index
  file <- resolved$file
  has_index <- !is.null(index)
  if (identical(mode, "indexed") && !has_index) {
    stop("`mode = \"indexed\"` requires a non-NULL `index`.", call. = FALSE)
  }
  if (has_index) {
    index <- validate_index(index)
  }

  files <- if (has_index) {
    if (is.null(file)) {
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
  } else {
    validate_input_files(file)
  }

  source_type <- if (has_index) {
    index$format
  } else {
    resolve_source_type_no_index(files, type)
  }

  use_streaming <- extraction_uses_streaming(
    mode,
    has_index,
    tuning$extract_mode
  )

  if (is.list(seq_idx)) {
    out <- vector("list", length(seq_idx))
    for (i in seq_along(seq_idx)) {
      out[[i]] <- extract_sequences_dnastringset(
        index = index,
        seq_idx = seq_idx[[i]],
        file = file,
        renumber = renumber,
        chunk_chars = chunk_chars,
        mode = mode,
        type = type
      )
    }
    return(out)
  }

  n_records <- if (has_index) {
    as.numeric(index$n_records)
  } else {
    roff <- streaming_record_offsets(NULL, files, source_type)
    as.numeric(roff[[length(roff)]])
  }
  if (length(seq_idx) < 1L) {
    return(Biostrings::DNAStringSet(character()))
  }
  seq_idx <- validate_seq_idx(seq_idx = seq_idx, n_records = n_records)

  if (use_streaming) {
    offsets <- streaming_record_offsets(index, files, source_type)
    # nolint nextline: object_usage_linter
    return(
      cpp_extract_sequences_dnastringset_streaming(
        files = files,
        source_type = source_type,
        ids_zero_based = as.numeric(seq_idx - 1),
        record_offsets_r = offsets,
        chunk_chars = as.numeric(chunk_chars),
        diagnostics = tuning$diagnostics,
        renumber_mode = renumber
      )
    )
  }

  live_index <- ensure_live_index_ptr(index)
  index <- live_index$index
  index_ptr <- live_index$ptr
  # nolint nextline: object_usage_linter
  cpp_extract_sequences_dnastringset(
    files = files,
    source_type = index$format,
    ids_zero_based = as.numeric(seq_idx - 1),
    index_ptr_sexp = index_ptr,
    chunk_chars = as.numeric(chunk_chars),
    max_bridge_gap = tuning$max_bridge_gap,
    max_region_records = tuning$max_region_records,
    extract_mode = tuning$extract_mode,
    diagnostics = tuning$diagnostics,
    renumber_mode = renumber
  )
}
