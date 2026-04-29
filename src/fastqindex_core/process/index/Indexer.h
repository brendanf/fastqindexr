/**
 * Adapted from FastqIndEx/src/process/index/Indexer.h (MIT).
 */

#ifndef FASTQINDEXR_INDEXER_H
#define FASTQINDEXR_INDEXER_H

#include "../base/IndexEntry.h"
#include "../base/ZLibBasedFASTQProcessorBaseClass.h"

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace fastqindex_core {

class Indexer : public ZLibBasedFASTQProcessorBaseClass {
 public:
  FileIndex createIndex(const std::string& gz_file_path);

 private:
  bool checkAndPrepareForNextConcatenatedPart(
    std::ifstream* in,
    z_stream* strm,
    std::vector<unsigned char>* window,
    std::vector<unsigned char>* input,
    bool* first_pass,
    bool* last_block_ended_with_newline,
    std::string* current_block,
    std::uint64_t total_bytes_in
  );

  void finalizeProcessingForCurrentBlock(
    const std::string& current_block,
    std::uint64_t block_offset,
    bool* first_pass,
    bool* last_block_ended_with_newline,
    int* cur_bits,
    std::uint64_t* line_count_for_next_index_entry,
    std::int64_t* block_id,
    z_stream* strm,
    std::vector<unsigned char>* dictionary_for_next_block,
    std::vector<IndexEntry>* entries
  );
};

// fastqindexr note (omitted upstream methods from Indexer):
// - fulfillsPremises()
// - createHeader()
// - checkAndPrepareForNextConcatenatedPart()
// - storeLinesOfCurrentBlockForDebugMode()
// - writeIndexEntryIfPossible() / storeDictionaryForEntry() split from file writer concerns
// - createIndexEntryFromBlockData() (logic in finalizeProcessingForCurrentBlock()).

}  // namespace fastqindex_core

#endif
