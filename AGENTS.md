# Agent / developer handoff: `fastqindexr`

This note is for future work in this repository (humans and coding agents). It captures intent, structure, and tooling quirks so you can stay aligned with how the project is meant to evolve.

## What this package is

- **R package** (MIT) that **indexes gzip-compressed FASTA or FASTQ** and **extracts records by position** through an **Rcpp** layer.
- The heavy lifting comes from a **vendored subset** of [FastqIndEx](https://github.com/DKFZ-ODCF/FastqIndEx) (MIT), adapted for in-memory use inside R—not a runtime dependency on upstream binaries.
- **Not** affiliated with the FastqIndEx maintainers; attribute upstream in `inst/LICENSE.note` and preserve copyright headers in vendored files.

## Design philosophy (do not fight these without an explicit product decision)

1. **Vendor upstream shape, adapt behavior.** C++ under `src/fastqindex_core/` keeps **the same directory layout and filenames** as FastqIndEx’s `src/` tree (`common/`, `process/base/`, `process/index/`, `process/extract/`). When refreshing from upstream, diff against that structure; avoid flattening or renaming “for convenience.”
2. **In-memory index, no `.fqi` on disk in the core path.** The R workflow builds and uses index data in process; the adapted core notes this where upstream assumed file-based index I/O.
3. **R-facing semantics:** record IDs in R are **1-based**; the C++ side uses 0-based offsets where appropriate (`extract_sequences` subtracts 1 before calling the bridge).
4. **Multiple files = one logical stream** for ID purposes: `create_index(c(f1, f2, ...))` assigns **global** record IDs in order; `file_record_offsets` in the index list defines which global IDs belong to which file.
5. **Output order = input order:** `extract_sequences(index, ids = c(3, 1, 3, ...))` returns rows in **that** order, including duplicates.
6. **Performance idea:** extraction batches requests into **dense regions** (see `kDenseGap` in `src/fastqindexr.cpp`) to reduce seek/decompression work, then reorders to the caller’s ID order.
7. **Known FASTA limitation (documented in Roxygen):** extraction assumes **one sequence line per FASTA record** (header + one line of sequence). Multi-line sequence blocks are not the target use case in the current implementation.
8. **Resilience:** if indexed random access fails for a region, the bridge can fall back to a **sequential** gzip parse for that file (see comments in `src/fastqindexr.cpp`).
9. **Single vendored tree only:** keep **one** FastqIndEx-derived source tree under `src/fastqindex_core/`. For new features, merge required upstream files into this tree (and reconcile same-name files there) instead of creating parallel vendored trees like `src/*_io` or `src/*_core2`.

## Code layout (where to look)

| Area | Role |
|------|------|
| `R/create_index.R`, `R/extract_sequences.R` | Public API, validation, S3 `print` for the index object |
| `R/imports.R` | `useDynLib` + `importFrom(Rcpp, evalCpp)` for NAMESPACE |
| `src/fastqindexr.cpp` | Rcpp exports, index registry, dense-region planning, parsing helpers |
| `src/fastqindex_core_bridge.cpp` | **Unity build:** `#include` of the vendored `.cpp` files **once** (avoids duplicate symbols) |
| `src/fastqindex_core/**` | Vendored FastqIndEx-style sources; `fastqindex_core` **namespace** isolates symbols |
| `src/Makevars` | `PKG_LIBS = -lz` (zlib) |
| `inst/LICENSE.note` | Upstream reference commit + list of vendored files + high-level change notes |
| `src/fastqindex_core/README.md` | Shorter note on layout and the bridge; keep in sync when the file set changes |
| `tests/testthat/` | Correctness tests; `tests/benchmarks/` for optional Biostrings comparisons |

Internal helpers `validate_input_files` / `validate_index` are **`@noRd`**: not part of the public API.

## C++ and Rcpp workflow

- After changing `[[Rcpp::export]]` signatures in `src/fastqindexr.cpp`, run **`Rcpp::compileAttributes()`** or **`pkgbuild::compile_dll()`** so `R/RcppExports.R` and `src/RcppExports.cpp` stay in sync.
- **Do not** compile the vendored `.cpp` files as separate translation units *in addition* to the bridge, or you will get duplicate symbol errors. The bridge is the single inclusion point.
- Vendored headers that use **quoted** includes (e.g. same-directory `#include "CommonStructsAndConstants.h"`) **rely** on the unity include order from the bridge; keep includes coherent.

## Documentation (Roxygen)

- **`DESCRIPTION`** sets `Roxygen: list(markdown = TRUE)`: write **markdown-style** Roxygen in `R/*.R` (backticks, `**bold**`, `[fun()]` links).
- Regenerate with **`roxygen2::roxygenize()`** (or `devtools::document()`). `NAMESPACE` is **roxygen-generated**; `R/imports.R` holds `useDynLib` / Rcpp imports.
- **`roxygen2` is in Suggests.** For markdown processing, a working **`commonmark`** install is expected in practice; without it, `roxygenize()` can fail or emit awkward Rd. If documentation looks wrong, check that dependency chain (historically `stringi` / `vctrs` in **renv** can also break loading if install was interrupted—see below).

## R environment: `renv`

- **`.Rprofile` sources `renv/activate.R`:** the project uses [renv](https://rstudio.github.io/renv/) for a project library.
- Agents may see **`renv::status()` “out of sync”** if `renv.lock` and installed packages differ—normal during dependency edits; run `renv::snapshot()` when the maintainer wants the lockfile updated.
- **Broken or half-installed packages** under the renv cache (missing `DESCRIPTION` in a package directory) can break *unrelated* loads (e.g. `roxygen2` failing because `stringi` is broken). Fix by reinstalling the broken package or cleaning the bad directory—not a problem with this repo’s R code per se.

## Continuous integration

- **`.github/workflows/r-cmd-check.yaml`:** `R CMD check` with `--no-manual` and CRAN-like args.
- **`.github/workflows/test-coverage.yaml`:** tests + **`covr::codecov()`**; repository must have Codecov (or similar) set up for badges to be meaningful.
- Pushes/PRs to **`main` or `master`**.

## README

- **User-facing** narrative lives in **`README.Rmd`**; `README.md` is generated—edit the `.Rmd` and knit if you change the public story.

## What to avoid in routine changes

- Reformatting or “cleaning up” vendored files beyond what a clear bug or upstream sync requires.
- **Drive-by** refactors unrelated to the task, or new markdown files unless the task asks for them.
- Renaming `src/fastqindex_core` tree paths to match personal taste—keep **upstream** alignment.
- Creating a second parallel vendored source tree for new upstream subsets; extend the existing `src/fastqindex_core/` tree instead.

## Quick verification commands (from the package root)

With renv active and dev dependencies available:

```r
pkgbuild::compile_dll()   # or devtools::load_all()
testthat::test_local()
roxygen2::roxygenize()
```

For a full sanity check: `R CMD build .` then `R CMD check` on the tarball.

---

*If you add or remove vendored files, update `inst/LICENSE.note` and `src/fastqindex_core/README.md`, and consider refreshing the **upstream commit** reference in `inst/LICENSE.note` to match the revision you actually compared or imported.*
