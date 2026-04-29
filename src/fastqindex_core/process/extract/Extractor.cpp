/**
 * Adapted from FastqIndEx/src/process/extract/Extractor.cpp (MIT).
 */

#include "Extractor.h"

#include <algorithm>
#ifdef FASTQINDEXR_TIMING
#include <chrono>
#endif
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace fastqindex_core {

// fastqindexr change (adapted from upstream Extractor::prepareForNextConcatenatedPartIfNecessary):
// Reinitialize raw inflate state after Z_STREAM_END so extraction can continue
// into concatenated gzip members when more lines are requested.
bool prepareForNextConcatenatedMember(
  std::ifstream* in,
  z_stream* strm,
  std::vector<unsigned char>* window,
  std::vector<unsigned char>* input,
  bool* first_pass,
  std::uint64_t total_bytes_in
) {
  if (in == nullptr || strm == nullptr || window == nullptr || input == nullptr || first_pass == nullptr) {
    return false;
  }
  // fastqindexr change: follow upstream offset rule (`+8+10`) before attempting
  // to continue in a concatenated gzip member.
  const std::uint64_t stream_end_position = total_bytes_in + 18;
  in->seekg(static_cast<std::streamoff>(stream_end_position), std::ios::beg);
  if (!in->good()) {
    return false;
  }

  inflateEnd(strm);
  std::memset(strm, 0, sizeof(z_stream));
  if (inflateInit2(strm, -15) != Z_OK) {
    return false;
  }
  std::fill(window->begin(), window->end(), 0);
  std::fill(input->begin(), input->end(), 0);
  strm->avail_out = WINDOW_SIZE;
  strm->next_out = window->data();
  unsigned char dict[WINDOW_SIZE]{0};
  if (inflateSetDictionary(strm, dict, WINDOW_SIZE) != Z_OK) {
    return false;
  }
  *first_pass = true;
  return true;
}

std::string trimLinePrefixLocal(const std::string& line, char prefix) {
  if (!line.empty() && line.front() == prefix) {
    return line.substr(1);
  }
  return line;
}

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
  std::uint64_t total_bytes_in = initial_offset;
  std::string incomplete_last_line;
  bool first_pass = true;
  std::vector<std::string> out;
  out.reserve(line_count);

  while (extracted_lines < line_count) {
    std::uint64_t total_bytes_in_dummy = 0;
    if (!readCompressedDataFromSource(in, &strm, input.data(), &total_bytes_in_dummy)) {
      break;
    }

    bool stream_ended = false;
    do {
      if (strm.avail_out == 0) {
        strm.avail_out = WINDOW_SIZE;
        strm.next_out = window.data();
      }

      std::uint64_t before_out = strm.avail_out;
      std::uint64_t before_in = strm.avail_in;
      int zlib_result = inflate(&strm, Z_NO_FLUSH);
      const std::uint64_t read_bytes = before_in - strm.avail_in;
      total_bytes_in += read_bytes;
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
        stream_ended = true;
        break;
      }
      if (zlib_result == Z_NEED_DICT || zlib_result == Z_DATA_ERROR || zlib_result == Z_MEM_ERROR) {
        inflateEnd(&strm);
        std::ostringstream oss;
        oss
          << "zlib extraction error (code=" << zlib_result
          << ", msg=" << (strm.msg == nullptr ? "none" : strm.msg)
          << ", start_line=" << starting_line
          << ", line_count=" << line_count
          << ", entry_start_line=" << entry->starting_line_in_entry
          << ", block_offset=" << entry->block_offset_in_raw_file
          << ", bits=" << static_cast<int>(entry->bits)
          << ", dict_compressed_size=" << static_cast<int>(entry->compressed_dictionary_size)
          << ").";
        throw std::runtime_error(oss.str());
      }
    } while (strm.avail_in != 0 && extracted_lines < line_count);

    if (stream_ended && extracted_lines < line_count) {
      if (!prepareForNextConcatenatedMember(
            &in,
            &strm,
            &window,
            &input,
            &first_pass,
            total_bytes_in
          )) {
        break;
      }
      total_bytes_in += 18;
    }
  }

  inflateEnd(&strm);
  return out;
}

void Extractor::extractSelected(
  const std::string& gz_file_path,
  const std::vector<IndexEntry>& entries,
  std::uint64_t starting_line,
  std::uint64_t line_count,
  int record_size,
  const std::unordered_set<std::uint64_t>& requested_local_ids,
  std::uint64_t region_start_record,
  const std::function<void(std::uint64_t, const std::vector<std::string>&)>& on_record
  #ifdef FASTQINDEXR_TIMING
  ,
  SelectiveExtractTimings* timings
  #endif
) {
#ifdef FASTQINDEXR_TIMING
  const auto setup_t0 = std::chrono::steady_clock::now();
#endif
  if (record_size <= 0) {
    throw std::runtime_error("Invalid record size for selective extraction.");
  }
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
  std::uint64_t processed_lines = 0;
  std::uint64_t total_bytes_in = initial_offset;
  std::string current_line;
  current_line.reserve(256);
  bool first_pass = true;
  std::uint64_t emitted_records = 0;
  std::size_t line_in_record = 0;
  std::vector<std::uint64_t> requested_sorted(requested_local_ids.begin(), requested_local_ids.end());
  std::sort(requested_sorted.begin(), requested_sorted.end());
  std::size_t requested_idx = 0;
  auto is_requested = [&](std::uint64_t local_id) -> bool {
    while (requested_idx < requested_sorted.size() && requested_sorted[requested_idx] < local_id) {
      requested_idx++;
    }
    return requested_idx < requested_sorted.size() && requested_sorted[requested_idx] == local_id;
  };
  bool capture_current_record = is_requested(region_start_record + emitted_records);
  std::vector<std::string> current_record;
  current_record.reserve(static_cast<std::size_t>(record_size));
  #ifdef FASTQINDEXR_TIMING
  if (timings != nullptr) {
#ifdef FASTQINDEXR_TIMING
    const auto setup_t1 = std::chrono::steady_clock::now();
    timings->time_seek_init_ms +=
      std::chrono::duration<double, std::milli>(setup_t1 - setup_t0).count();
#endif
  }
  #endif

  while (processed_lines < line_count) {
    std::uint64_t total_bytes_in_dummy = 0;
    if (!readCompressedDataFromSource(in, &strm, input.data(), &total_bytes_in_dummy)) {
      break;
    }

    bool stream_ended = false;
    do {
      if (strm.avail_out == 0) {
        strm.avail_out = WINDOW_SIZE;
        strm.next_out = window.data();
      }

      std::uint64_t before_out = strm.avail_out;
      std::uint64_t before_in = strm.avail_in;
#ifdef FASTQINDEXR_TIMING
      const auto inflate_t0 = std::chrono::steady_clock::now();
#endif
      int zlib_result = inflate(&strm, Z_NO_FLUSH);
      #ifdef FASTQINDEXR_TIMING
      if (timings != nullptr) {
#ifdef FASTQINDEXR_TIMING
        const auto inflate_t1 = std::chrono::steady_clock::now();
        timings->time_inflate_ms +=
          std::chrono::duration<double, std::milli>(inflate_t1 - inflate_t0).count();
#endif
      }
      #endif
      const std::uint64_t read_bytes = before_in - strm.avail_in;
      total_bytes_in += read_bytes;
      std::uint64_t written = before_out - strm.avail_out;
      const unsigned char* chunk_start = window.data() + (WINDOW_SIZE - before_out);
#ifdef FASTQINDEXR_TIMING
      const auto parse_t0 = std::chrono::steady_clock::now();
      const auto split_t0 = std::chrono::steady_clock::now();
#endif
      auto consume_line = [&](const std::string& line) {
        // fastqindexr change: preserve first-pass behavior from the prior
        // split-based implementation while parsing in-place to avoid per-chunk
        // line-vector allocations.
        if (first_pass && entry->offset_to_next_line_start > 0) {
          first_pass = false;
          if (capture_current_record) {
            current_line.clear();
          }
          return;
        }
        first_pass = false;

        if (skip > 0) {
          skip--;
          if (capture_current_record) {
            current_line.clear();
          }
          return;
        }
        if (processed_lines >= line_count) {
          if (capture_current_record) {
            current_line.clear();
          }
          return;
        }
        if (capture_current_record) {
#ifdef FASTQINDEXR_TIMING
          const auto materialize_t0 = std::chrono::steady_clock::now();
#endif
          current_record.push_back(line);
          #ifdef FASTQINDEXR_TIMING
          if (timings != nullptr) {
#ifdef FASTQINDEXR_TIMING
            const auto materialize_t1 = std::chrono::steady_clock::now();
            timings->time_materialize_ms +=
              std::chrono::duration<double, std::milli>(materialize_t1 - materialize_t0).count();
#endif
          }
          #endif
        }
        processed_lines++;
        line_in_record++;
        if (line_in_record == static_cast<std::size_t>(record_size)) {
          const std::uint64_t local_id = region_start_record + emitted_records;
          if (capture_current_record) {
#ifdef FASTQINDEXR_TIMING
            const auto callback_t0 = std::chrono::steady_clock::now();
#endif
            on_record(local_id, current_record);
            #ifdef FASTQINDEXR_TIMING
            if (timings != nullptr) {
#ifdef FASTQINDEXR_TIMING
              const auto callback_t1 = std::chrono::steady_clock::now();
              timings->time_callback_ms +=
                std::chrono::duration<double, std::milli>(callback_t1 - callback_t0).count();
#endif
            }
            #endif
            current_record.clear();
          }
          emitted_records++;
          line_in_record = 0;
          capture_current_record = is_requested(region_start_record + emitted_records);
        }
      };

      std::size_t line_start = 0;
      for (std::size_t i = 0; i < static_cast<std::size_t>(written); ++i) {
        if (chunk_start[i] != '\n') {
          continue;
        }
        std::size_t seg_end = i;
        if (seg_end > line_start && chunk_start[seg_end - 1] == '\r') {
          seg_end--;
        }
        if (capture_current_record) {
          if (seg_end > line_start) {
            const char* seg_ptr = reinterpret_cast<const char*>(chunk_start + line_start);
            const std::size_t seg_len = seg_end - line_start;
            if (current_line.empty()) {
              current_line.assign(seg_ptr, seg_len);
            } else {
              current_line.append(seg_ptr, seg_len);
            }
          }
          consume_line(current_line);
          current_line.clear();
        } else {
          consume_line(std::string{});
        }
        line_start = i + 1;
      }
      if (line_start < static_cast<std::size_t>(written) && capture_current_record) {
        const char* seg_ptr = reinterpret_cast<const char*>(chunk_start + line_start);
        const std::size_t seg_len = static_cast<std::size_t>(written) - line_start;
        if (current_line.empty()) {
          current_line.assign(seg_ptr, seg_len);
        } else {
          current_line.append(seg_ptr, seg_len);
        }
      }
      #ifdef FASTQINDEXR_TIMING
      if (timings != nullptr) {
#ifdef FASTQINDEXR_TIMING
        const auto split_t1 = std::chrono::steady_clock::now();
        timings->time_line_split_ms +=
          std::chrono::duration<double, std::milli>(split_t1 - split_t0).count();
        const auto parse_t1 = std::chrono::steady_clock::now();
        timings->time_parse_ms +=
          std::chrono::duration<double, std::milli>(parse_t1 - parse_t0).count();
#endif
      }
      #endif

      if (zlib_result == Z_STREAM_END) {
        stream_ended = true;
        break;
      }
      if (zlib_result == Z_NEED_DICT || zlib_result == Z_DATA_ERROR || zlib_result == Z_MEM_ERROR) {
        inflateEnd(&strm);
        std::ostringstream oss;
        oss
          << "zlib extraction error (code=" << zlib_result
          << ", msg=" << (strm.msg == nullptr ? "none" : strm.msg)
          << ", start_line=" << starting_line
          << ", line_count=" << line_count
          << ", entry_start_line=" << entry->starting_line_in_entry
          << ", block_offset=" << entry->block_offset_in_raw_file
          << ", bits=" << static_cast<int>(entry->bits)
          << ", dict_compressed_size=" << static_cast<int>(entry->compressed_dictionary_size)
          << ").";
        throw std::runtime_error(oss.str());
      }
    } while (strm.avail_in != 0 && processed_lines < line_count);

    if (stream_ended && processed_lines < line_count) {
      if (!prepareForNextConcatenatedMember(
            &in,
            &strm,
            &window,
            &input,
            &first_pass,
            total_bytes_in
          )) {
        break;
      }
      total_bytes_in += 18;
    }
  }

  inflateEnd(&strm);
}

void Extractor::extractSelectedRecords(
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
  ,
  SelectiveExtractTimings* timings
  #endif
) {
  const bool is_fastq = (type == "fastq");
#ifdef FASTQINDEXR_TIMING
  const auto setup_t0 = std::chrono::steady_clock::now();
#endif
  if (record_size <= 0) {
    throw std::runtime_error("Invalid record size for selective extraction.");
  }
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
  std::uint64_t processed_lines = 0;
  std::uint64_t total_bytes_in = initial_offset;
  std::string current_line;
  current_line.reserve(256);
  bool first_pass = true;
  std::uint64_t emitted_records = 0;
  std::size_t line_in_record = 0;
  bool capture_current_record =
    requested_local_ids.find(region_start_record + emitted_records) != requested_local_ids.end();
  ExtractedRecord current_record;
  #ifdef FASTQINDEXR_TIMING
  if (timings != nullptr) {
#ifdef FASTQINDEXR_TIMING
    const auto setup_t1 = std::chrono::steady_clock::now();
    timings->time_seek_init_ms +=
      std::chrono::duration<double, std::milli>(setup_t1 - setup_t0).count();
#endif
  }
  #endif

  while (processed_lines < line_count) {
    std::uint64_t total_bytes_in_dummy = 0;
    if (!readCompressedDataFromSource(in, &strm, input.data(), &total_bytes_in_dummy)) {
      break;
    }

    bool stream_ended = false;
    do {
      if (strm.avail_out == 0) {
        strm.avail_out = WINDOW_SIZE;
        strm.next_out = window.data();
      }

      std::uint64_t before_out = strm.avail_out;
      std::uint64_t before_in = strm.avail_in;
#ifdef FASTQINDEXR_TIMING
      const auto inflate_t0 = std::chrono::steady_clock::now();
#endif
      int zlib_result = inflate(&strm, Z_NO_FLUSH);
      #ifdef FASTQINDEXR_TIMING
      if (timings != nullptr) {
#ifdef FASTQINDEXR_TIMING
        const auto inflate_t1 = std::chrono::steady_clock::now();
        timings->time_inflate_ms +=
          std::chrono::duration<double, std::milli>(inflate_t1 - inflate_t0).count();
#endif
      }
      #endif
      const std::uint64_t read_bytes = before_in - strm.avail_in;
      total_bytes_in += read_bytes;
      std::uint64_t written = before_out - strm.avail_out;
      const unsigned char* chunk_start = window.data() + (WINDOW_SIZE - before_out);
#ifdef FASTQINDEXR_TIMING
      const auto parse_t0 = std::chrono::steady_clock::now();
      const auto split_t0 = std::chrono::steady_clock::now();
#endif
      for (std::uint64_t i = 0; i < written; ++i) {
        const unsigned char ch = chunk_start[i];
        if (ch == '\r') {
          continue;
        }
        if (ch != '\n') {
          if (capture_current_record) {
            current_line.push_back(static_cast<char>(ch));
          }
          continue;
        }
        if (first_pass && entry->offset_to_next_line_start > 0) {
          first_pass = false;
          if (capture_current_record) {
            current_line.clear();
          }
          continue;
        }
        first_pass = false;
        if (skip > 0) {
          skip--;
          if (capture_current_record) {
            current_line.clear();
          }
          continue;
        }
        if (processed_lines >= line_count) {
          if (capture_current_record) {
            current_line.clear();
          }
          break;
        }
        if (capture_current_record) {
#ifdef FASTQINDEXR_TIMING
          const auto materialize_t0 = std::chrono::steady_clock::now();
#endif
          if (is_fastq) {
            if (line_in_record == 0) {
              current_record.seq_id = trimLinePrefixLocal(current_line, '@');
            } else if (line_in_record == 1) {
              current_record.seq = current_line;
            } else if (line_in_record == 3) {
              current_record.qual = include_qual ? current_line : "";
            }
          } else {
            if (line_in_record == 0) {
              current_record.seq_id = trimLinePrefixLocal(current_line, '>');
            } else if (line_in_record == 1) {
              current_record.seq = current_line;
            }
          }
          current_line.clear();
          #ifdef FASTQINDEXR_TIMING
          if (timings != nullptr) {
#ifdef FASTQINDEXR_TIMING
            const auto materialize_t1 = std::chrono::steady_clock::now();
            timings->time_materialize_ms +=
              std::chrono::duration<double, std::milli>(materialize_t1 - materialize_t0).count();
#endif
          }
          #endif
        }
        processed_lines++;
        line_in_record++;
        if (line_in_record == static_cast<std::size_t>(record_size)) {
          const std::uint64_t local_id = region_start_record + emitted_records;
          if (capture_current_record) {
#ifdef FASTQINDEXR_TIMING
            const auto callback_t0 = std::chrono::steady_clock::now();
#endif
            on_record(local_id, current_record);
            #ifdef FASTQINDEXR_TIMING
            if (timings != nullptr) {
#ifdef FASTQINDEXR_TIMING
              const auto callback_t1 = std::chrono::steady_clock::now();
              timings->time_callback_ms +=
                std::chrono::duration<double, std::milli>(callback_t1 - callback_t0).count();
#endif
            }
            #endif
            current_record = ExtractedRecord{};
          }
          emitted_records++;
          line_in_record = 0;
          capture_current_record =
            requested_local_ids.find(region_start_record + emitted_records) != requested_local_ids.end();
        }
      }
      #ifdef FASTQINDEXR_TIMING
      if (timings != nullptr) {
#ifdef FASTQINDEXR_TIMING
        const auto split_t1 = std::chrono::steady_clock::now();
        timings->time_line_split_ms +=
          std::chrono::duration<double, std::milli>(split_t1 - split_t0).count();
        const auto parse_t1 = std::chrono::steady_clock::now();
        timings->time_parse_ms +=
          std::chrono::duration<double, std::milli>(parse_t1 - parse_t0).count();
#endif
      }
      #endif

      if (zlib_result == Z_STREAM_END) {
        stream_ended = true;
        break;
      }
      if (zlib_result == Z_NEED_DICT || zlib_result == Z_DATA_ERROR || zlib_result == Z_MEM_ERROR) {
        inflateEnd(&strm);
        throw std::runtime_error("zlib extraction error during selective record parsing.");
      }
    } while (strm.avail_in != 0 && processed_lines < line_count);

    if (stream_ended && processed_lines < line_count) {
      if (!prepareForNextConcatenatedMember(
            &in,
            &strm,
            &window,
            &input,
            &first_pass,
            total_bytes_in
          )) {
        break;
      }
      total_bytes_in += 18;
    }
  }
  inflateEnd(&strm);
}

void Extractor::extractRecords(
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
  ,
  SelectiveExtractTimings* timings
  #endif
) {
  std::unordered_set<std::uint64_t> all_requested;
  const std::uint64_t n_records = line_count / static_cast<std::uint64_t>(record_size);
  all_requested.reserve(static_cast<std::size_t>(n_records));
  for (std::uint64_t i = 0; i < n_records; ++i) {
    all_requested.insert(region_start_record + i);
  }
  extractSelectedRecords(
    gz_file_path,
    entries,
    starting_line,
    line_count,
    record_size,
    all_requested,
    region_start_record,
    type,
    include_qual,
    on_record
    #ifdef FASTQINDEXR_TIMING
    , timings
    #endif
  );
}

}  // namespace fastqindex_core
