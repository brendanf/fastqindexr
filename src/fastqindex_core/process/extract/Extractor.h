/**
 * Adapted from FastqIndEx/src/process/extract/Extractor.h (MIT).
 */

#ifndef FASTQINDEXR_EXTRACTOR_H
#define FASTQINDEXR_EXTRACTOR_H

#include "../base/IndexEntry.h"
#include "../base/ZLibBasedFASTQProcessorBaseClass.h"

#include <cstdint>
#include <string>
#include <vector>

namespace fastqindex_core {

class Extractor : public ZLibBasedFASTQProcessorBaseClass {
 public:
  std::vector<std::string> extract(
    const std::string& gz_file_path,
    const std::vector<IndexEntry>& entries,
    std::uint64_t starting_line,
    std::uint64_t line_count
  );

 private:
  const IndexEntry* findIndexEntryForExtraction(
    const std::vector<IndexEntry>& entries,
    std::uint64_t starting_line
  ) const;
};

// fastqindexr note (omitted upstream methods from Extractor):
// - fulfillsPremises(), calculateStartingLineAndLineCount()
// - openFastqAndPrepareZStream(), setDictionaryForZStream()
// - prepareForNextConcatenatedPartIfNecessary()
// - processDecompressedChunkOfData(), storeOrOutputLine(), storeLinesOfCurrentBlockForDebugMode()
// fastqindexr exposes direct line-window extraction for R-managed record planning.

}  // namespace fastqindex_core

#endif
