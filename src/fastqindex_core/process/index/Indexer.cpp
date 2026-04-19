/**
 * Adapted from FastqIndEx/src/process/index/Indexer.cpp (MIT).
 */

#include "Indexer.h"

#include <cstring>
#include <fstream>
#include <stdexcept>

namespace fastqindex_core {

void Indexer::finalizeProcessingForCurrentBlock(
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
) {
  if (*first_pass) {
    *first_pass = false;
    return;
  }

  *block_id += 1;
  std::vector<std::string> lines = splitStr(current_block);
  std::uint32_t number_of_lines = static_cast<std::uint32_t>(lines.size());
  bool current_block_ended_with_newline =
    current_block.empty() ? *last_block_ended_with_newline : current_block.back() == '\n';

  bool has_any_line_breaks = current_block.find('\n') != std::string::npos;
  std::uint32_t offset_first_line = 0;
  if (!*last_block_ended_with_newline) {
    if (number_of_lines > 0) {
      number_of_lines -= 1;
    }
    if (has_any_line_breaks && !lines.empty()) {
      offset_first_line = static_cast<std::uint32_t>(lines.front().size() + 1);
    }
  }

  IndexEntry entry;
  entry.block_index = static_cast<std::uint64_t>(*block_id);
  entry.block_offset_in_raw_file = block_offset;
  entry.starting_line_in_entry = *line_count_for_next_index_entry;
  entry.offset_to_next_line_start = offset_first_line;
  entry.bits = static_cast<unsigned char>(*cur_bits);
  entry.dictionary = *dictionary_for_next_block;

  // fastqindexr change: we keep dictionary uncompressed in-memory instead of writing .fqi.
  std::vector<unsigned char> next_dictionary(WINDOW_SIZE, 0);
  unsigned int copied_bytes = 0;
  inflateGetDictionary(strm, next_dictionary.data(), &copied_bytes);
  *dictionary_for_next_block = std::move(next_dictionary);

  entries->push_back(std::move(entry));
  *line_count_for_next_index_entry += number_of_lines;
  *last_block_ended_with_newline = current_block_ended_with_newline;
  *cur_bits = strm->data_type & 7;
}

FileIndex Indexer::createIndex(const std::string& gz_file_path) {
  // fastqindexr change: follows upstream block-boundary indexing flow but stores
  // entries in-memory (FileIndex) instead of using IndexWriter/.fqi output.
  std::ifstream in(gz_file_path, std::ios::binary);
  if (!in.good()) {
    throw std::runtime_error("Could not open compressed file: " + gz_file_path);
  }

  z_stream strm;
  if (!initializeZStream(&strm, 47)) {
    throw std::runtime_error("Could not initialize zlib stream for indexing.");
  }

  std::vector<unsigned char> input(CHUNK_SIZE, 0);
  std::vector<unsigned char> window(WINDOW_SIZE, 0);
  std::vector<unsigned char> dictionary_for_next_block(WINDOW_SIZE, 0);
  strm.avail_out = WINDOW_SIZE;
  strm.next_out = window.data();

  bool keep_processing = true;
  bool first_pass = true;
  bool last_block_ended_with_newline = true;
  int cur_bits = 0;
  std::int64_t block_id = -1;
  std::uint64_t total_bytes_in = 0;
  std::uint64_t offset = 0;
  std::uint64_t line_count_for_next_index_entry = 0;
  std::string current_block;
  std::vector<IndexEntry> entries;

  while (keep_processing) {
    std::uint64_t block_offset = offset;
    if (!readCompressedDataFromSource(in, &strm, input.data(), &total_bytes_in)) {
      break;
    }
    offset = total_bytes_in;

    do {
      if (strm.avail_out == 0) {
        strm.avail_out = WINDOW_SIZE;
        strm.next_out = window.data();
      }
      std::uint64_t before_out = strm.avail_out;
      int zlib_result = inflate(&strm, Z_BLOCK);
      std::uint64_t written = before_out - strm.avail_out;
      current_block.append(reinterpret_cast<char*>(window.data() + (WINDOW_SIZE - before_out)), written);

      if (zlib_result == Z_NEED_DICT || zlib_result == Z_DATA_ERROR || zlib_result == Z_MEM_ERROR) {
        inflateEnd(&strm);
        throw std::runtime_error("zlib indexing error.");
      }

      if (zlib_result == Z_STREAM_END) {
        finalizeProcessingForCurrentBlock(
          current_block,
          block_offset,
          &first_pass,
          &last_block_ended_with_newline,
          &cur_bits,
          &line_count_for_next_index_entry,
          &block_id,
          &strm,
          &dictionary_for_next_block,
          &entries
        );
        keep_processing = false;
        break;
      }

      if (checkStreamForBlockEnd(strm)) {
        finalizeProcessingForCurrentBlock(
          current_block,
          block_offset,
          &first_pass,
          &last_block_ended_with_newline,
          &cur_bits,
          &line_count_for_next_index_entry,
          &block_id,
          &strm,
          &dictionary_for_next_block,
          &entries
        );
        current_block.clear();
        block_offset = total_bytes_in - strm.avail_in;
      }
    } while (strm.avail_in != 0);
  }

  inflateEnd(&strm);
  FileIndex out;
  out.entries = std::move(entries);
  out.total_lines = line_count_for_next_index_entry;
  return out;
}

}  // namespace fastqindex_core
