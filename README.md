
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastqindexr

<!-- badges: start -->

[![R-CMD-check](https://github.com/brendanf/fastqindexr/actions/workflows/r-cmd-check.yaml/badge.svg)](https://github.com/brendanf/fastqindexr/actions/workflows/r-cmd-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/brendanf/fastqindexr/graph/badge.svg)](https://app.codecov.io/gh/brendanf/fastqindexr)
<!-- badges: end -->

`fastqindexr` builds an in-memory index over one or more gzipped
FASTA/FASTQ files, then extracts records by ID without scanning from the
beginning each time. It wraps an Rcpp bridge around adapted open-source
[FastqIndEx](https://github.com/DKFZ-ODCF/FastqIndEx) C++ core
components. It is not written or maintained by the authors of
FastqIndEx.

The package is useful when you want to:

- repeatedly sample or subset records from large `.gz` inputs
- treat several files as one logical concatenated record stream
- keep extracted output in R as a data frame (`seq_id`, `seq`, and
  `qual` for FASTQ)

## Installation

You can install the development version of fastqindexr from GitHub with:

``` r
# install.packages("pak")
pak::pak("brendanf/fastqindexr")
```

## Example

Build an index once, then request records in any order (including
duplicates):

``` r
library(fastqindexr)

# Make a tiny gzipped FASTQ for demonstration.
path <- tempfile(fileext = ".fastq.gz")
con <- gzfile(path, "wt")
writeLines(
  c(
    "@read1", "ACGT", "+", "!!!!",
    "@read2", "TTAA", "+", "####",
    "@read3", "GGCC", "+", "$$$$"
  ),
  con
)
close(con)

# Index the file (type can also be "auto" or "fasta").
idx <- create_index(path, type = "fastq")
idx
#> <fastqindexr_index>
#>   format: fastq 
#>   files: 1 
#>   records: 3

# Extract by 1-based record ID, preserving request order.
extract_sequences(idx, seq_idx = c(3, 1, 1))
#>   seq_id  seq qual
#> 1  read3 GGCC $$$$
#> 2  read1 ACGT !!!!
#> 3  read1 ACGT !!!!

unlink(path)
```

## Benchmarking

The following benchmarks compares `fastqindexr` against the `Biostrings`

``` r
# Install Biostrings if not already installed.
# install.packages("BiocManager")
# BiocManager::install("Biostrings")

set.seed(1)

tmp <- tempfile(fileext = ".fa.gz")
make_benchmark_fasta(tmp, n = 10000, width = 120)
ids <- sample.int(10000, 2000, replace = TRUE)

# fastqindexr index creation
system.time(
  idx <- create_index(tmp, type = "fasta")
)
#>    user  system elapsed 
#>   0.009   0.001   0.008

# Biostrings index creation
system.time(
  bi_index <- Biostrings::fasta.index(tmp, seqtype = "DNA")
)
#>    user  system elapsed 
#>   1.387   0.073   1.485

# fastqindexr indexed extraction
system.time({
  res_fastqindexr <- extract_sequences(idx, ids)
})
#>    user  system elapsed 
#>    0.01    0.00    0.01

# Biostrings indexed extraction
system.time({
  selected <- bi_index[ids, , drop = FALSE]
  res_biostrings <- Biostrings::readDNAStringSet(selected)
})
#>    user  system elapsed 
#>   3.635   0.052   3.693

# verify that sequences and named are identical
all.equal(
  res_fastqindexr$seq,
  as.character(res_biostrings),
  check.attributes = FALSE
)
#> [1] TRUE

all.equal(
  res_fastqindexr$seq_id,
  names(res_biostrings)
)
#> [1] TRUE

unlink(tmp)
```

## API at a glance

- `create_index(files, type = c("auto", "fasta", "fastq"))`
  - accepts one or many existing gzipped files
  - returns a `fastqindexr_index` object
  - this object can be saved to disk with `saveRDS()` and loaded in a
    different session with `readRDS()` (or serialized/deserialized in
    other ways, such as with the
    [`qs2`](https://cran.r-project.org/package=qs2) package)
- `extract_sequences(index, seq_idx, file = NULL)`
  - `seq_idx` are 1-based positive integer record IDs
  - returns rows in the same order as `seq_idx`
  - duplicate IDs are allowed and duplicated in output
  - for FASTQ, output includes `seq_id`, `seq`, and `qual`; for FASTA,
    it includes `seq_id` and `seq`
- `extract_sequences_to_file(index, seq_idx, file = NULL, outfile, ...)`
  - streams extracted records directly to file (plain or `.gz`)
  - preserves input order and duplicate IDs exactly
  - supports `type = "auto"`, `"fasta"`, or `"fastq"` output
- `make_benchmark_fasta(path, n = 5000, width = 80)`
  - writes synthetic gzipped FASTA for repeatable benchmark setup

## Important format assumptions

- FASTA support currently assumes one sequence line per record (header
  line + one sequence line)
- files are treated as a logical concatenated stream in the order
  supplied to `create_index()`
- if files move after indexing, you can provide replacement paths via
  the `file` argument in `extract_sequences()`
