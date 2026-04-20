#' Write a synthetic gzipped FASTA file for benchmarking
#'
#' Creates a FASTA file where each record has a deterministic ID
#' (`seq00001`, `seq00002`, ...) and a randomly sampled DNA sequence.
#'
#' @param path Output file path, typically ending in `.fa.gz` or `.fasta.gz`.
#' @param n Number of FASTA records to write.
#' @param width Sequence width (bases) per record.
#' @param alphabet Character vector of symbols used to generate sequence lines.
#'
#' @return Invisibly returns `path`.
#'
#' @examples
#' path <- tempfile(fileext = ".fa.gz")
#' make_benchmark_fasta(path, n = 3, width = 10)
#' file.exists(path)
#' unlink(path)
#'
#' @export
make_benchmark_fasta <- function(
  path,
  n = 5000,
  width = 80,
  alphabet = c("A", "C", "G", "T")
) {
  path <- as.character(path)
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single non-empty string.", call. = FALSE)
  }

  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1) {
    stop("`n` must be a single positive number.", call. = FALSE)
  }
  n <- as.integer(n)

  if (!is.numeric(width) || length(width) != 1L || is.na(width) || width < 1) {
    stop("`width` must be a single positive number.", call. = FALSE)
  }
  width <- as.integer(width)

  alphabet <- as.character(alphabet)
  if (length(alphabet) < 1L || anyNA(alphabet) || any(!nzchar(alphabet))) {
    stop(
      "`alphabet` must contain at least one non-empty symbol.",
      call. = FALSE
    )
  }

  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)

  for (i in seq_len(n)) {
    seq <- paste0(sample(alphabet, width, replace = TRUE), collapse = "")
    writeLines(c(sprintf(">seq%05d", i), seq), con = con)
  }

  invisible(path)
}
