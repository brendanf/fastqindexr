# fastqindexr

## Development

- Add `extract_sequences_to_file()` to stream extracted records directly to
  plain-text or gzip-compressed output files.
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