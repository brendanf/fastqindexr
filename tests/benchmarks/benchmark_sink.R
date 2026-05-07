# Shared helpers for benchmark scripts: tee console output to a small git-trackable
# snapshot under tests/benchmarks/results/ for regression diffing.

bench_results_dir <- function() {
  d <- file.path("tests", "benchmarks", "results")
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  d
}

bench_sink_start <- function(script_name) {
  path <- file.path(bench_results_dir(), paste0(script_name, ".txt"))
  git_sha <- tryCatch(
    trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)[1]),
    error = function(...) "unknown"
  )
  writeLines(
    c(
      paste("# fastqindexr benchmark snapshot:", script_name),
      paste("# git:", git_sha),
      paste("# R:", R.version.string),
      "---",
      ""
    ),
    path
  )
  sink(file(path, open = "at"), split = TRUE, type = "output")
  invisible(path)
}

bench_sink_stop <- function() {
  while (sink.number() > 0L) {
    sink()
  }
}

bench_log_header <- function(label) {
  cat(sprintf("\n[benchmark]\n%s\n", label))
}

bench_log_timing <- function(label, timing) {
  cat(sprintf("case: %s\n", label))
  cat(sprintf("elapsed_s: %.3f\n", timing[["elapsed"]]))
}
