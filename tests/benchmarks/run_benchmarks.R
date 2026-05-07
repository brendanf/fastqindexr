for (file in list.files("tests/benchmarks", pattern = "benchmark.R$")) {
  cat(sprintf("Running %s\n", file))
  source(file)
}