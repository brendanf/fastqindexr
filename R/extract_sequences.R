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
      !(extract_mode %in% c("indexed", "sequential_only"))
  ) {
    stop(
      paste0(
        "Option `fastqindexr.extract_mode` must be one of ",
        "`\"indexed\"` or `\"sequential_only\"`."
      ),
      call. = FALSE
    )
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
  tuning <- resolve_region_merge_tuning()
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

  n_records <- as.numeric(index$n_records)
  seq_idx <- validate_seq_idx(seq_idx = seq_idx, n_records = n_records)

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
    include_qual = include_qual,
    max_bridge_gap = tuning$max_bridge_gap,
    max_region_records = tuning$max_region_records,
    extract_mode = tuning$extract_mode,
    diagnostics = tuning$diagnostics
  )
  diag <- attr(result, "fastqindexr_diagnostics", exact = TRUE)

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
  index,
  seq_idx,
  file = NULL,
  outfile,
  type = c("auto", "fasta", "fastq"),
  append = FALSE,
  compress = endsWith(tolower(outfile), ".gz")
) {
  tuning <- resolve_region_merge_tuning()
  resolved <- resolve_extract_index(index, file)
  index <- resolved$index
  file <- resolved$file
  index <- validate_index(index)
  live_index <- ensure_live_index_ptr(index)
  index <- live_index$index
  index_ptr <- live_index$ptr

  if (!is.logical(append) || length(append) != 1L || is.na(append)) {
    stop("`append` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(outfile) || length(outfile) < 1L || any(!nzchar(outfile))) {
    stop("`outfile` must be a non-empty character vector.", call. = FALSE)
  }
  is_partition_mode <- is.list(seq_idx) || length(outfile) > 1L

  n_records <- as.numeric(index$n_records)

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
      resolved_type <- if (type == "auto") index$format else type
      if (index$format == "fasta" && resolved_type == "fastq") {
        stop("Cannot emit FASTQ output from FASTA index input.", call. = FALSE)
      }
      include_qual <- identical(index$format, "fastq") &&
        identical(resolved_type, "fastq")
      union_ids <- sort(all_ids)
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
  index,
  seq_idx,
  file = NULL,
  renumber = c("none", "zero_based", "one_based"),
  chunk_chars = 1e7
) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package `Biostrings` is required.", call. = FALSE)
  }
  renumber <- match.arg(renumber)
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
  index <- validate_index(resolved$index)
  live_index <- ensure_live_index_ptr(index)
  index <- live_index$index
  index_ptr <- live_index$ptr
  n_records <- as.numeric(index$n_records)
  if (length(seq_idx) < 1L) {
    return(Biostrings::DNAStringSet(character()))
  }
  seq_idx <- validate_seq_idx(seq_idx = seq_idx, n_records = n_records)
  files <- if (is.null(resolved$file)) {
    as.character(index$files)
  } else {
    files <- validate_input_files(resolved$file)
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
