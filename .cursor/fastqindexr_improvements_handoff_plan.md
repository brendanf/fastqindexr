# fastqindexr Improvements Handoff Plan

## Scope

Implement optional `fastqindexr` enhancements identified during
`optimotu.pipeline` migration planning, without disrupting current exported API
or vendored-core alignment.

Primary targets:

- Add an extraction path that writes directly to file/connection to reduce peak
  memory and copies for large requests.
- Add extraction output modes to avoid unnecessary data.frame overhead when only
  sequence payloads are needed.
- Add an optimized fast path for sorted/near-contiguous IDs.

## Non-negotiable constraints from project guidance

- Keep the vendored FastqIndEx tree structure under `src/fastqindex_core/`.
- Do not introduce a second vendored tree.
- Preserve current semantics:
  - 1-based record IDs at R boundary
  - multi-file logical concatenation behavior
  - output order identical to requested ID order, including duplicates
- Keep known FASTA one-line-per-record limitation unless explicitly expanding
  scope.
- If vendored files change, update:
  - `inst/LICENSE.note`
  - `src/fastqindex_core/README.md`

## Workstream A: Streaming extraction API

### Goal

Avoid materializing full extracted result frames in R when caller wants output
written to file.

### Proposed public API

Add new function:

- `extract_sequences_to_file(index, seq_idx, file = NULL, outfile, type = c("auto", "fasta", "fastq"), append = FALSE, compress = endsWith(outfile, ".gz"))`

Notes:

- Keep existing `extract_sequences()` unchanged.
- `index` / `file` semantics should mirror `extract_sequences()` exactly.
- `type` can usually be inferred from index metadata; keep argument mostly for
  explicit output control and future-proofing.

### Implementation outline

1. Add R wrapper in `R/extract_sequences.R` with validation shared with
   `extract_sequences()`.
2. Add C++ bridge entrypoint in `src/fastqindexr.cpp` that:
   - reuses index pointer restoration machinery
   - extracts records in requested order
   - writes FASTA/FASTQ lines directly to output stream
3. Support gzipped output path via either:
   - zlib stream in C++, or
   - safe staged write plus gzip pipe, whichever is simpler and robust.
4. Preserve duplicate IDs and ordering exactly.

### Tests

- New tests in `tests/testthat/test-extract-sequences-to-file.R`:
  - FASTA and FASTQ outputs
  - duplicate IDs and arbitrary order
  - append and overwrite behavior
  - object index vs `.fqi` index path
  - multi-file index path vector

## Workstream B: Return-mode variants for in-memory extraction

### Goal

Reduce overhead when consumers do not need full data.frame payloads.

### Proposed API approach

Option 1 (preferred for minimal API growth):

- extend `extract_sequences()` with `return = c("data.frame", "list", "seq")`

Semantics:

- `data.frame` (default): current behavior
- `list`: list with fields (`seq_id`, `seq`, and `qual` for FASTQ)
- `seq`: character vector of sequences only, with names as `seq_id`

### Implementation outline

1. Keep current default untouched to avoid downstream breakage.
2. Implement post-processing in R first (lower risk), then optionally move
   specialized paths into C++ if profiling shows clear wins.
3. Ensure output ordering and duplication invariants remain unchanged.

### Tests

- Add coverage in `tests/testthat/test-extract-sequences.R`:
  - all return modes for FASTA and FASTQ
  - edge case `length(seq_idx) == 0`
  - duplicate IDs and name retention for `return = "seq"`

## Workstream C: Sorted-ID fast path

### Goal

Improve throughput for callers that provide sorted and locally dense IDs.

### Implementation outline

1. In `src/fastqindexr.cpp`, add detection:
   - already sorted non-decreasing IDs
   - contiguous/dense runs above threshold
2. Reuse existing dense-region planning (`kDenseGap`) but bypass unnecessary
   reorder bookkeeping when input is already sorted and no duplicates require
   remapping.
3. Keep correctness-first fallback to current generic path.

### Benchmarks and safety checks

- Add benchmark script under `tests/benchmarks/` comparing:
  - random IDs
  - sorted dense IDs
  - sorted sparse IDs
- Validate no regression for random-ID path.

## Documentation and build updates

1. Update roxygen docs in `R/extract_sequences.R`.
2. Regenerate docs and namespace with roxygen2:
   - do not hand-edit `NAMESPACE` or `.Rd`.
3. If new exports are added, ensure `README.Rmd` API section includes them.
4. If changing Rcpp exports, run `Rcpp::compileAttributes()`.

## Validation checklist for handoff agent

1. `pkgbuild::compile_dll()`
2. `testthat::test_local()`
3. `roxygen2::roxygenize()`
4. Re-run focused extraction tests after doc regeneration.
5. Optional: `R CMD build` + `R CMD check` tarball.

## Suggested execution order

1. Streaming extraction API (highest impact for memory).
2. Return-mode variants (low-risk usability/perf win).
3. Sorted-ID fast path (profiling-guided optimization).

## Out-of-scope for this handoff

- Expanding FASTA support to multi-line sequence records.
- Structural reorganization of vendored core directories.
- Replacing the current C++ bridge architecture.
