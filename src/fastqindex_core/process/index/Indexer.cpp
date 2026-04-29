/**
 * Adapted from FastqIndEx/src/process/index/Indexer.cpp (MIT).
 */

#include "Indexer.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <stdexcept>

namespace fastqindex_core {

bool Indexer::checkAndPrepareForNextConcatenatedPart(
  std::ifstream* in,
  z_stream* strm,
  std::vector<unsigned char>* window,
  std::vector<unsigned char>* input,
  bool* first_pass,
  bool* last_block_ended_with_newline,
  std::string* current_block,
  std::uint64_t total_bytes_in
) {
  if (in == nullptr || strm == nullptr || window == nullptr || input == nullptr ||
      first_pass == nullptr || last_block_ended_with_newline == nullptr || current_block == nullptr) {
    return false;
  }

  // fastqindexr change: copied from upstream Indexer::checkAndPrepareForNextConcatenatedPart()
  // and adapted for std::ifstream/std::vector storage.
  in->seekg(static_cast<std::streamoff>(total_bytes_in), std::ios::beg);
  if (!in->good()) {
    return false;
  }

  inflateEnd(strm);
  if (!initializeZStream(strm, 47)) {
    return false;
  }
  std::fill(window->begin(), window->end(), 0);
  std::fill(input->begin(), input->end(), 0);
  strm->avail_out = WINDOW_SIZE;
  strm->next_out = window->data();
  strm->avail_in = CHUNK_SIZE;
  *first_pass = true;
  *last_block_ended_with_newline = true;
  current_block->clear();
  return true;
}

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
  entry.compressed_dictionary_size = 0;
  entry.dictionary = *dictionary_for_next_block;

  // fastqindexr change: we keep dictionary uncompressed in-memory instead of writing .fqi.
  // fastqindexr change: mirror upstream behavior by reusing the same backing
  // buffer for `inflateGetDictionary()`; this preserves any untouched prefix
  // bytes exactly as in FastqIndEx's fixed-size array implementation.
  unsigned int copied_bytes = 0;
  inflateGetDictionary(strm, dictionary_for_next_block->data(), &copied_bytes);

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
  int zlib_result = Z_OK;
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
    do {
      if (!in.good()) {
        keep_processing = false;
        break;
      }

      std::uint64_t total_bytes_in_dummy = 0;
      if (!readCompressedDataFromSource(in, &strm, input.data(), &total_bytes_in_dummy)) {
        keep_processing = false;
        break;
      }

      do {
        if (strm.avail_out == 0) {
          strm.avail_out = WINDOW_SIZE;
          strm.next_out = window.data();
        }
        const std::uint64_t before_in = strm.avail_in;
        std::uint64_t before_out = strm.avail_out;
        zlib_result = inflate(&strm, Z_BLOCK);
        const std::uint64_t read_bytes = before_in - strm.avail_in;
        total_bytes_in += read_bytes;
        std::uint64_t written = before_out - strm.avail_out;
        current_block.append(reinterpret_cast<char*>(window.data() + (WINDOW_SIZE - before_out)), written);

        if (zlib_result == Z_NEED_DICT || zlib_result == Z_DATA_ERROR || zlib_result == Z_MEM_ERROR) {
          inflateEnd(&strm);
          throw std::runtime_error("zlib indexing error.");
        }

        if (zlib_result == Z_STREAM_END) {
          finalizeProcessingForCurrentBlock(
            current_block,
            offset,
            &first_pass,
            &last_block_ended_with_newline,
            &cur_bits,
            &line_count_for_next_index_entry,
            &block_id,
            &strm,
            &dictionary_for_next_block,
            &entries
          );
          break;
        }

        if (checkStreamForBlockEnd(strm)) {
          finalizeProcessingForCurrentBlock(
            current_block,
            offset,
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
          offset = total_bytes_in;
        }

        // fastqindexr change: mirror upstream control flow and clear first-pass
        // after each inflate iteration.
        first_pass = false;
      } while (strm.avail_in != 0);

      keep_processing = zlib_result != Z_STREAM_END && in.good();
    } while (keep_processing);

    // fastqindexr change: upstream parity for concatenated gzip continuation.
    if (!keep_processing) {
      keep_processing = checkAndPrepareForNextConcatenatedPart(
        &in,
        &strm,
        &window,
        &input,
        &first_pass,
        &last_block_ended_with_newline,
        &current_block,
        total_bytes_in
      );
      if (keep_processing) {
        offset = total_bytes_in;
      }
    }
  }

  inflateEnd(&strm);
  FileIndex out;
  out.entries = std::move(entries);
  out.total_lines = line_count_for_next_index_entry;
  return out;
}

}  // namespace fastqindex_core
