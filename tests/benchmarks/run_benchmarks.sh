# Run all benchmarks in the tests/benchmarks directory.
# From package root:
# tests/benchmarks/run_benchmarks.sh

for file in tests/benchmarks/*benchmark.R; do
  echo "Running $file"
  Rscript "$file"
done
