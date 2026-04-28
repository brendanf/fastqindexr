/**
 * Adapted from FastqIndEx/src/process/extract/Extractor.cpp (MIT).
 */

#include "Extractor.h"

#include <algorithm>
#include <fstream>
#include <stdexcept>

namespace fastqindex_core {

const IndexEntry* Extractor::findIndexEntryForExtraction(
  const std::vector<IndexEntry>& entries,
  std::uint64_t starting_line
) const {
  if (entries.empty()) {
    return nullptr;
  }
  const IndexEntry* latest = &entries.front();
  for (const auto& entry : entries) {
    if (entry.starting_line_in_entry > starting_line) {
      break;
    }
    latest = &entry;
  }
  return latest;
}

std::vector<std::string> Extractor::extract(
  const std::string& gz_file_path,
  const std::vector<IndexEntry>& entries,
  std::uint64_t starting_line,
  std::uint64_t line_count
) {
  // fastqindexr change: preserves upstream random-seek inflatePrime/inflateSetDictionary
  // logic, but returns decompressed lines in-memory instead of writing to a Sink.
  const IndexEntry* entry = findIndexEntryForExtraction(entries, starting_line);
  if (entry == nullptr) {
    throw std::runtime_error("Index contains no entries.");
  }

  std::ifstream in(gz_file_path, std::ios::binary);
  if (!in.good()) {
    throw std::runtime_error("Could not open compressed file: " + gz_file_path);
  }

  z_stream strm;
  if (!initializeZStream(&strm, -15)) {
    throw std::runtime_error("Could not initialize zlib raw stream for extraction.");
  }

  std::uint64_t initial_offset = entry->block_offset_in_raw_file;
  unsigned char start_bits = entry->bits;
  if (start_bits > 0) {
    if (initial_offset == 0) {
      inflateEnd(&strm);
      throw std::runtime_error("Invalid index entry offset/bits combination.");
    }
    initial_offset -= 1;
  }
  in.seekg(static_cast<std::streamoff>(initial_offset), std::ios::beg);
  if (!in.good()) {
    inflateEnd(&strm);
    throw std::runtime_error("Could not seek to indexed offset.");
  }

  if (start_bits > 0) {
    int next = in.get();
    if (next == EOF) {
      inflateEnd(&strm);
      throw std::runtime_error("Could not read primer byte for inflatePrime.");
    }
    if (inflatePrime(&strm, start_bits, next >> (8 - start_bits)) != Z_OK) {
      inflateEnd(&strm);
      throw std::runtime_error("inflatePrime failed.");
    }
  }

  std::vector<unsigned char> dictionary_for_inflate;
  const unsigned char* dictionary_ptr = nullptr;
  if (entry->compressed_dictionary_size > 0) {
    uLongf dest_len = WINDOW_SIZE;
    dictionary_for_inflate.assign(WINDOW_SIZE, 0);
    uLong source_len = entry->compressed_dictionary_size;
    if (entry->dictionary.size() < source_len ||
        uncompress2(
          dictionary_for_inflate.data(),
          &dest_len,
          entry->dictionary.data(),
          &source_len
        ) != Z_OK ||
        dest_len != WINDOW_SIZE) {
      inflateEnd(&strm);
      throw std::runtime_error("Failed to decompress index dictionary.");
    }
    dictionary_ptr = dictionary_for_inflate.data();
  } else if (entry->dictionary.size() == WINDOW_SIZE) {
    dictionary_ptr = entry->dictionary.data();
  } else {
    inflateEnd(&strm);
    throw std::runtime_error("Unexpected dictionary size in index entry.");
  }

  if (inflateSetDictionary(&strm, dictionary_ptr, WINDOW_SIZE) != Z_OK) {
    inflateEnd(&strm);
    throw std::runtime_error("inflateSetDictionary failed.");
  }

  std::vector<unsigned char> input(CHUNK_SIZE, 0);
  std::vector<unsigned char> window(WINDOW_SIZE, 0);
  strm.avail_out = WINDOW_SIZE;
  strm.next_out = window.data();

  std::uint64_t skip = starting_line - entry->starting_line_in_entry;
  std::uint64_t extracted_lines = 0;
  std::string incomplete_last_line;
  bool first_pass = true;
  std::vector<std::string> out;
  out.reserve(line_count);

  while (extracted_lines < line_count) {
    std::uint64_t total_bytes_in_dummy = 0;
    if (!readCompressedDataFromSource(in, &strm, input.data(), &total_bytes_in_dummy)) {
      break;
    }

    do {
      if (strm.avail_out == 0) {
        strm.avail_out = WINDOW_SIZE;
        strm.next_out = window.data();
      }

      std::uint64_t before_out = strm.avail_out;
      int zlib_result = inflate(&strm, Z_NO_FLUSH);
      std::uint64_t written = before_out - strm.avail_out;
      std::string chunk(reinterpret_cast<char*>(window.data() + (WINDOW_SIZE - before_out)), written);

      // fastqindexr change: chunk handling is adapted from processDecompressedChunkOfData()
      // and kept local to avoid depending on Sink/error-accumulator subsystems.
      std::vector<std::string> split = splitStr(chunk);
      if (first_pass && entry->offset_to_next_line_start > 0 && !split.empty()) {
        split.erase(split.begin());
      }
      first_pass = false;

      std::string cur_incomplete;
      if (!chunk.empty() && chunk.back() != '\n' && !split.empty()) {
        cur_incomplete = split.back();
        split.pop_back();
      }

      if (!split.empty()) {
        if (skip >= split.size()) {
          skip -= split.size();
        } else {
          std::size_t i = static_cast<std::size_t>(skip);
          if (i < split.size()) {
            split[i] = incomplete_last_line + split[i];
            for (; i < split.size() && extracted_lines < line_count; ++i) {
              out.push_back(split[i]);
              extracted_lines++;
            }
          }
          skip = 0;
        }
      }
      incomplete_last_line = cur_incomplete;

      if (zlib_result == Z_STREAM_END) {
        break;
      }
      if (zlib_result == Z_NEED_DICT || zlib_result == Z_DATA_ERROR || zlib_result == Z_MEM_ERROR) {
        inflateEnd(&strm);
        throw std::runtime_error("zlib extraction error.");
      }
    } while (strm.avail_in != 0 && extracted_lines < line_count);
  }

  inflateEnd(&strm);
  return out;
}

}  // namespace fastqindex_core
