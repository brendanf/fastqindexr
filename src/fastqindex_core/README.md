This directory holds vendored/adapted FastqIndEx core C++ sources.

Source project: https://github.com/DKFZ-ODCF/FastqIndEx (MIT License)

The same tree now contains both in-memory indexing/extraction adaptations and
the vendored `.fqi` reader subset (no separate second source tree).

**Layout:** Paths and filenames mirror the upstream `src/` tree: `common/`, `process/base/`,
`process/index/`, `process/extract/`, etc., rooted here as `src/fastqindex_core/...` (not
`src/src/...`). The R package adds `fastqindex_core_bridge.cpp` beside this folder; it
includes the vendored `.cpp` files once (unity compilation) so translation units link
exactly once. Types used from R live in the `fastqindex_core` namespace to avoid global
symbol clashes.

When copying or refreshing upstream files:
1. Preserve original copyright/license headers.
2. Record upstream commit/tag and local modifications in `inst/LICENSE.note`.

Current adapted files:
- `common/CommonStructsAndConstants.h`
- `common/ErrorAccumulator.h`
- `common/ErrorAccumulator.cpp`
- `common/IOHelper.h`
- `common/IOHelper.cpp`
- `process/base/IndexEntry.h`
- `process/base/BaseIndexEntry.h`
- `process/base/IndexEntryV1.h`
- `process/base/IndexHeader.h`
- `process/base/ZLibBasedFASTQProcessorBaseClass.h`
- `process/base/ZLibBasedFASTQProcessorBaseClass.cpp`
- `process/io/IOBase.h`
- `process/io/Source.h`
- `process/io/FileSource.h`
- `process/io/FileSource.cpp`
- `process/io/locks/LockHandler.h`
- `process/io/locks/FileLockHandler.h`
- `process/io/locks/FileLockHandler.cpp`
- `process/index/Indexer.h`
- `process/index/Indexer.cpp`
- `process/extract/IndexReader.h`
- `process/extract/IndexReader.cpp`
- `process/extract/Extractor.h`
- `process/extract/Extractor.cpp`

Each adapted file contains inline notes for:
- `fastqindexr change` blocks in modified functions used by this package.
- omitted upstream functions that are intentionally not included in this adapted subset.
