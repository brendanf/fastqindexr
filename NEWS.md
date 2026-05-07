# fastqindexr 0.1.0

- Add extraction `mode = c("auto", "indexed", "sequential")` across
  `extract_sequences()`, `extract_sequences_to_file()`, and
  `extract_sequences_dnastringset()`. Keep
  `fastqindexr.extract_mode` as a transition option and deprecate
  `"sequential_only"` in favor of `"sequential"`.
- Support sequential extraction without a live native pointer: with
  `index = NULL`, extraction streams directly from source file(s)
  (`file`, `type`), and with `mode = "sequential"` extraction can run
  from serialized index metadata after `saveRDS()` / `readRDS()`.
- Extend index creation and indexed extraction to uncompressed FASTA/FASTQ
  alongside gzip inputs.
- Refactor index internals into gzip/plain subtype hierarchies and require
  homogeneous compression in each index (mixed gzip/plain inputs now error).
- Replace serialized payload schema with the new typed format
  (`schema_version` v4); older serialized payloads are no longer restored
  and must be rebuilt.
- Add `partition_seq_idx()` and extend partition-aware extraction:
  `extract_sequences_to_file()` supports partitioned writes, and
  `extract_sequences()` / `extract_sequences_dnastringset()` accept
  list-valued `seq_idx` and return partition-aligned results.
- Add `renumber = c("none", "zero_based", "one_based")` to
  `extract_sequences()` and `extract_sequences_to_file()`, and align
  empty-request behavior across in-memory outputs and partitioned files.
- Improve extraction performance and robustness: major sparse-extraction
  speedups, chunked memory behavior for `extract_sequences()` and
  `extract_sequences_dnastringset()`, region-merge tuning via
  `fastqindexr.max_bridge_gap` / `fastqindexr.max_region_records`,
  and fixes for some FastqIndEx `.fqi` reads and repeated small
  extractions with reused indexes.
- Support **multi-line FASTA sequences** for both gzip and plain inputs
  across indexed extraction, sequential extraction, and file output.
- Add `collapse_sequence_lines` to `extract_sequences_to_file()` for
  FASTA output, allowing sequence lines to be emitted either with their
  original line breaks (default) or collapsed to a single line.
- Improve plain-file FASTA performance via offset-bounded indexed reads
  and lower-overhead sequential parsing for skipped records.
- Extend benchmark helpers with `make_benchmark_fastq()` alongside
  `make_benchmark_fasta()` for synthetic FASTQ generation.
- Update serialized index payload handling to support schema version 5
  (while still restoring schema version 4 payloads).

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
