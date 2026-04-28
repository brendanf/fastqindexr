# fastqindexr (development version)

- Fix a bug reading certain .fqi indexes written by FastqIndEx.
- Fix a bug when reusing the same index for many small extractions.
- Add `partition_seq_idx()` for stable contiguous or round-robin partitioning
  of requested sequence IDs while preserving order and duplicates.
- Extend `extract_sequences_to_file()` with partitioned mode by passing a list
  to `seq_idx` and a matching vector to `outfile`. Supports scalar/vector
  `compress` and preserves strict-increasing fast-path behavior per partition.
- Add `extract_sequences_dnastringset()` with chunked construction of
  `Biostrings::DNAStringSet` for improved large-request memory behavior.
- Add similar chunked extraction to `extract_sequences()` for reduced memory usage in large extracts.

## 0.0.2

- Faster extraction when `seq_idx` is **strictly increasing** (so 0-based global
  IDs are sorted with no duplicates): skip per-file `sort` / `unique` on local
  record indices. Unsorted or duplicate `seq_idx` uses the previous path.
- Add `extract_sequences_to_file()` to stream extracted records directly to
  plain-text or gzip-compressed output files. For **strictly increasing**
  `seq_idx` with `append = FALSE`, file output is written as each record is
  read (no in-memory map of all unique records). `append = TRUE` uses a
  buffered path like in-memory extraction.
- Extend `extract_sequences()` with `return = c("data.frame", "list", "seq")`
  for list or named-sequence output without the full `data.frame` when not
  needed.
- Skip reading and returning quality when not needed: FASTA-indexed
  extraction omits a `qual` field; `return = "seq"` for FASTQ omits
  per-base quality; `extract_sequences_to_file(..., type = "fasta")` with a
  FASTQ index does not read quality from disk.


## 0.0.1

Initial release. The package builds an in-memory index over one or more gzipped
FASTA or FASTQ files and extracts records by 1-based ID without a full pass
from the start of the file each time. The implementation uses an Rcpp bridge
around adapted [FastqIndEx](https://github.com/DKFZ-ODCF/FastqIndEx) C++
components (this package is not maintained by the FastqIndEx authors). Several
files are treated as one logical concatenated record stream. Main functions:
`create_index()`, `extract_sequences()` (data frame with `seq_id` and `seq`, and
`qual` for FASTQ), `read_fqi_index()` to use FastqIndEx CLI `.fqi` indexes, and
`make_benchmark_fasta()` for synthetic gzipped FASTA in benchmarks. Indexes can
be saved and reloaded (for example with `saveRDS()` / `readRDS()`).
