/**
 * Adapted from FastqIndEx/src/process/extract/Extractor.h (MIT).
 */

#ifndef FASTQINDEXR_EXTRACTOR_H
#define FASTQINDEXR_EXTRACTOR_H

#include "../base/IndexEntry.h"
#include "../base/ZLibBasedFASTQProcessorBaseClass.h"

#include <cstdint>
#include <functional>
#include <string>
#include <unordered_set>
#include <vector>

namespace fastqindex_core {

#ifdef FASTQINDEXR_TIMING
struct SelectiveExtractTimings {
  double time_seek_init_ms{0.0};
  double time_inflate_ms{0.0};
  double time_parse_ms{0.0};
  double time_line_split_ms{0.0};
  double time_materialize_ms{0.0};
  double time_callback_ms{0.0};
};
#endif

class Extractor : public ZLibBasedFASTQProcessorBaseClass {
 public:
  std::vector<std::string> extract(
    const std::string& gz_file_path,
    const std::vector<IndexEntry>& entries,
    std::uint64_t starting_line,
    std::uint64_t line_count
  );
  void extractSelected(
    const std::string& gz_file_path,
    const std::vector<IndexEntry>& entries,
    std::uint64_t starting_line,
    std::uint64_t line_count,
    int record_size,
    const std::unordered_set<std::uint64_t>& requested_local_ids,
    std::uint64_t region_start_record,
    const std::function<void(std::uint64_t, const std::vector<std::string>&)>& on_record
    #ifdef FASTQINDEXR_TIMING
    , SelectiveExtractTimings* timings = nullptr
    #endif
  );
  void extractSelectedRecords(
    const std::string& gz_file_path,
    const std::vector<IndexEntry>& entries,
    std::uint64_t starting_line,
    std::uint64_t line_count,
    int record_size,
    const std::unordered_set<std::uint64_t>& requested_local_ids,
    std::uint64_t region_start_record,
    const std::string& type,
    bool include_qual,
    const std::function<void(std::uint64_t, const ExtractedRecord&)>& on_record
    #ifdef FASTQINDEXR_TIMING
    , SelectiveExtractTimings* timings = nullptr
    #endif
  );
  void extractRecords(
    const std::string& gz_file_path,
    const std::vector<IndexEntry>& entries,
    std::uint64_t starting_line,
    std::uint64_t line_count,
    int record_size,
    std::uint64_t region_start_record,
    const std::string& type,
    bool include_qual,
    const std::function<void(std::uint64_t, const ExtractedRecord&)>& on_record
    #ifdef FASTQINDEXR_TIMING
    , SelectiveExtractTimings* timings = nullptr
    #endif
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
