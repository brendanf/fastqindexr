/**
 * Adapted from FastqIndEx/src/process/base/IndexEntry*.h (MIT).
 */

#ifndef FASTQINDEXR_INDEXENTRY_H
#define FASTQINDEXR_INDEXENTRY_H

#include <cstdint>
#include <string>
#include <vector>

namespace fastqindex_core {

struct IndexEntry {
  std::uint64_t block_index{0};
  std::uint64_t block_offset_in_raw_file{0};
  std::uint64_t starting_line_in_entry{0};
  std::uint32_t offset_to_next_line_start{0};
  unsigned char bits{0};
  std::vector<unsigned char> dictionary;
};

struct FileIndex {
  std::vector<IndexEntry> entries;
  std::uint64_t total_lines{0};
};

struct ExtractedRecord {
  std::string seq_id;
  std::string seq;
  std::string qual;
};

}  // namespace fastqindex_core

#endif
