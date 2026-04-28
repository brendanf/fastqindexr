# fastqindexr Handoff: Indexed External-Tool Helper Support

## Context

`optimotu.pipeline` is being refactored so external-tool wrappers can accept:

- a `fastqindexr_index` object (or `.fqi` path),
- optional `seq_idx`,
- optional `nthreads` (for batch partitioning),

and then own the full extract/split/process lifecycle internally.

The goal for `fastqindexr` is to provide a small set of reusable primitives that
reduce orchestration code in downstream wrappers while preserving current
correctness guarantees (ordering, duplicates, and multi-file logical indexing).

## Requested features

### 1) Stable partition helper for `seq_idx`

Add an exported helper that partitions requested sequence IDs into N batches.

Suggested signature:

```r
partition_seq_idx <- function(
  seq_idx,
  n_parts,
  strategy = c("contiguous", "round_robin"),
  drop_empty = TRUE
)
```

Expected behavior:

- accepts positive integer-like `seq_idx`; preserves original values.
- `strategy = "contiguous"`: split input vector into near-equal contiguous
  chunks, preserving in-chunk order.
- `strategy = "round_robin"`: assign elements by position modulo `n_parts`,
  preserving relative order within each partition.
- supports duplicates and out-of-order values.
- returns `list(integer())` partitions (length `n_parts` if `drop_empty = FALSE`,
  otherwise only non-empty).

Why:

- Downstream helpers currently reimplement batching logic ad hoc.
- A shared helper improves readability and testability.

### 2) Multi-output extraction in one call

Add an API to extract multiple `seq_idx` partitions directly to multiple output
files, reducing R-level loops and repeated setup overhead.

Suggested signature:

```r
extract_sequences_partitions_to_files <- function(
  index,
  seq_idx_partitions,
  file = NULL,
  outfiles,
  type = c("auto", "fasta", "fastq"),
  append = FALSE,
  compress = endsWith(tolower(outfiles), ".gz")
)
```

Expected behavior:

- `seq_idx_partitions` is a list of numeric vectors (each validated like
  `extract_sequences_to_file()`).
- `length(outfiles)` must match `length(seq_idx_partitions)`.
- each partition writes only its own IDs, preserving order and duplicates.
- returns normalized output paths invisibly.
- empty partition:
  - if `append = FALSE`, create/truncate corresponding output file as empty.
  - if `append = TRUE`, keep existing file untouched.
- supports index object and `.fqi` path mode as existing extract APIs do.

Why:

- `optimotu.pipeline` helpers often need "extract then run one process per
  partition"; this API makes stage 1+2 concise and consistent.

### 3) Direct DNAStringSet return helper (optional but preferred)

Add a convenience wrapper to return `DNAStringSet` directly from indexed inputs.

Suggested signature:

```r
extract_sequences_dnastringset <- function(
  index,
  seq_idx,
  file = NULL,
  renumber = c("none", "zero_based", "one_based")
)
```

Expected behavior:

- extracts sequence and IDs from FASTA/FASTQ index and returns
  `Biostrings::DNAStringSet`.
- names preserve original IDs when `renumber = "none"`.
- renumbering replaces names deterministically:
  - `"zero_based"` => `"0"`, `"1"`, ...
  - `"one_based"` => `"1"`, `"2"`, ...
- preserves request order and duplicates.

Why:

- downstream wrappers frequently convert extraction outputs to DNAStringSet.

## Correctness contracts (must hold)

For all new APIs:

- request order is preserved exactly.
- duplicate IDs are preserved.
- behavior is correct for multi-file logical indices.
- bounds validation errors match current strictness.
- empty requests return typed empty outputs (not `NULL`).

## Benchmark scenarios

Include at least these benchmark fixtures in package-level benchmarking or
repeatable examples:

1) Single file, ordered IDs:
   - `seq_idx = 1:200000`
2) Single file, random IDs with duplicates:
   - `seq_idx` random sample with replacement
3) Two-file logical index:
   - request spans file boundary and includes duplicates near boundary
4) Partitioned extraction:
   - N-way split and multi-output write for N in `{1, 2, 8, 32}`

Capture:

- wall time,
- peak memory (if feasible),
- output equivalence against current extract APIs.

## Acceptance criteria for downstream (`optimotu.pipeline`)

`optimotu.pipeline` migration depends on:

1) New helper API availability with docs + tests.
2) Deterministic behavior for duplicate/out-of-order `seq_idx`.
3) Multi-file index support in all new APIs.
4) No regression in existing `extract_sequences()` /
   `extract_sequences_to_file()` behavior.

## Notes for implementation scope

- Prefer composing existing validation internals where possible.
- Keep new features additive; avoid breaking old signatures.
- Roxygen docs should clearly state ordering/duplicate guarantees.

