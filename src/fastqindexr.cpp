#include <Rcpp.h>
#include <zlib.h>

#include "fastqindex_core/process/extract/Extractor.h"
#include "fastqindex_core/process/extract/IndexReader.h"
#include "fastqindex_core/process/index/Indexer.h"
#include "fastqindex_core/process/io/FileSource.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef WINDOW_SIZE
#undef WINDOW_SIZE
#endif

namespace {

constexpr int kGzBufferSize = 8192;
constexpr std::uint64_t kDenseGap = 64;

struct ParsedRecord {
  std::string id;
  std::string seq;
  std::string qual;
};

struct IndexedFile {
  fastqindex_core::FileIndex index;
  std::uint64_t record_count{0};
};

struct IndexBundle {
  std::string format;
  int record_size{0};
  std::vector<std::string> files;
  std::vector<IndexedFile> indexed_files;
  std::vector<std::uint64_t> record_offsets;
  std::uint64_t n_records{0};
};

class OutputWriter {
 public:
  virtual ~OutputWriter() = default;
  virtual void write(const std::string& chunk) = 0;
};

class PlainOutputWriter : public OutputWriter {
 public:
  PlainOutputWriter(const std::string& path, bool append) {
    std::ios_base::openmode mode = std::ios::binary | std::ios::out;
    mode |= append ? std::ios::app : std::ios::trunc;
    stream_.open(path, mode);
    if (!stream_.good()) {
      throw std::runtime_error("Could not open output file for writing: " + path);
    }
  }

  void write(const std::string& chunk) override {
    stream_.write(chunk.data(), static_cast<std::streamsize>(chunk.size()));
    if (!stream_.good()) {
      throw std::runtime_error("Write failed while streaming extraction output.");
    }
  }

 private:
  std::ofstream stream_;
};

class GzipOutputWriter : public OutputWriter {
 public:
  GzipOutputWriter(const std::string& path, bool append) {
    const char* mode = append ? "ab" : "wb";
    stream_ = gzopen(path.c_str(), mode);
    if (stream_ == nullptr) {
      throw std::runtime_error("Could not open gz output stream: " + path);
    }
  }

  ~GzipOutputWriter() override {
    if (stream_ != nullptr) {
      gzclose(stream_);
    }
  }

  void write(const std::string& chunk) override {
    const int n = static_cast<int>(chunk.size());
    if (gzwrite(stream_, chunk.data(), n) != n) {
      int zerr = Z_OK;
      const char* msg = gzerror(stream_, &zerr);
      std::string detail = (msg == nullptr) ? "unknown zlib error" : msg;
      throw std::runtime_error("gzwrite failed while streaming extraction output: " +
                               detail);
    }
  }

 private:
  gzFile stream_{nullptr};
};

using SelectedRecordMap = std::unordered_map<std::uint64_t, ParsedRecord>;

std::uint64_t asUInt64(double value, const std::string& field) {
  if (!std::isfinite(value) || value < 0) {
    throw std::runtime_error("Malformed serialized index field: " + field);
  }
  return static_cast<std::uint64_t>(value);
}

Rcpp::List serializeIndexedFile(const IndexedFile& indexed_file) {
  const auto& entries = indexed_file.index.entries;
  const R_xlen_t n = static_cast<R_xlen_t>(entries.size());

  Rcpp::NumericVector block_index(n);
  Rcpp::NumericVector block_offset(n);
  Rcpp::NumericVector starting_line(n);
  Rcpp::NumericVector offset_next_line(n);
  Rcpp::IntegerVector compressed_dictionary_size(n);
  Rcpp::IntegerVector bits(n);
  Rcpp::RawVector dictionary_blob(
    static_cast<R_xlen_t>(entries.size() * fastqindex_core::WINDOW_SIZE)
  );

  for (R_xlen_t i = 0; i < n; ++i) {
    const auto& entry = entries[static_cast<size_t>(i)];
    if (entry.dictionary.size() != fastqindex_core::WINDOW_SIZE) {
      throw std::runtime_error("Unexpected dictionary size in index entry.");
    }
    block_index[i] = static_cast<double>(entry.block_index);
    block_offset[i] = static_cast<double>(entry.block_offset_in_raw_file);
    starting_line[i] = static_cast<double>(entry.starting_line_in_entry);
    offset_next_line[i] = static_cast<double>(entry.offset_to_next_line_start);
    compressed_dictionary_size[i] = static_cast<int>(entry.compressed_dictionary_size);
    bits[i] = static_cast<int>(entry.bits);

    const size_t dictionary_offset = static_cast<size_t>(i) * fastqindex_core::WINDOW_SIZE;
    std::copy(
      entry.dictionary.begin(),
      entry.dictionary.end(),
      dictionary_blob.begin() + static_cast<R_xlen_t>(dictionary_offset)
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("block_index") = block_index,
    Rcpp::Named("block_offset_in_raw_file") = block_offset,
    Rcpp::Named("starting_line_in_entry") = starting_line,
    Rcpp::Named("offset_to_next_line_start") = offset_next_line,
    Rcpp::Named("compressed_dictionary_size") = compressed_dictionary_size,
    Rcpp::Named("bits") = bits,
    Rcpp::Named("dictionary_blob") = dictionary_blob,
    Rcpp::Named("total_lines") = static_cast<double>(indexed_file.index.total_lines),
    Rcpp::Named("record_count") = static_cast<double>(indexed_file.record_count)
  );
}

Rcpp::List serializeIndexBundle(const IndexBundle& bundle) {
  Rcpp::List indexed_files(bundle.indexed_files.size());
  for (R_xlen_t i = 0; i < indexed_files.size(); ++i) {
    indexed_files[i] = serializeIndexedFile(bundle.indexed_files[static_cast<size_t>(i)]);
  }

  Rcpp::NumericVector record_offsets(bundle.record_offsets.size());
  for (R_xlen_t i = 0; i < record_offsets.size(); ++i) {
    record_offsets[i] = static_cast<double>(bundle.record_offsets[static_cast<size_t>(i)]);
  }

  return Rcpp::List::create(
    Rcpp::Named("schema_version") = 2L,
    Rcpp::Named("format") = bundle.format,
    Rcpp::Named("record_size") = bundle.record_size,
    Rcpp::Named("files") = bundle.files,
    Rcpp::Named("record_offsets") = record_offsets,
    Rcpp::Named("n_records") = static_cast<double>(bundle.n_records),
    Rcpp::Named("indexed_files") = indexed_files
  );
}

IndexedFile deserializeIndexedFile(const Rcpp::List& payload, int schema_version) {
  Rcpp::NumericVector block_index = payload["block_index"];
  Rcpp::NumericVector block_offset = payload["block_offset_in_raw_file"];
  Rcpp::NumericVector starting_line = payload["starting_line_in_entry"];
  Rcpp::NumericVector offset_next_line = payload["offset_to_next_line_start"];
  Rcpp::IntegerVector compressed_dictionary_size;
  if (schema_version >= 2 && payload.containsElementNamed("compressed_dictionary_size")) {
    compressed_dictionary_size = payload["compressed_dictionary_size"];
  }
  Rcpp::IntegerVector bits = payload["bits"];
  Rcpp::RawVector dictionary_blob = payload["dictionary_blob"];

  const R_xlen_t n = block_offset.size();
  if (block_index.size() != n || starting_line.size() != n || offset_next_line.size() != n || bits.size() != n) {
    throw std::runtime_error("Malformed serialized indexed file payload.");
  }
  if (compressed_dictionary_size.size() > 0 && compressed_dictionary_size.size() != n) {
    throw std::runtime_error("Malformed compressed dictionary size payload.");
  }
  if (dictionary_blob.size() != static_cast<R_xlen_t>(n * fastqindex_core::WINDOW_SIZE)) {
    throw std::runtime_error("Malformed dictionary blob in serialized payload.");
  }

  IndexedFile indexed_file;
  indexed_file.index.entries.reserve(static_cast<size_t>(n));
  for (R_xlen_t i = 0; i < n; ++i) {
    fastqindex_core::IndexEntry entry;
    entry.block_index = asUInt64(block_index[i], "block_index");
    entry.block_offset_in_raw_file = asUInt64(block_offset[i], "block_offset_in_raw_file");
    entry.starting_line_in_entry = asUInt64(starting_line[i], "starting_line_in_entry");
    entry.offset_to_next_line_start = static_cast<std::uint32_t>(asUInt64(offset_next_line[i], "offset_to_next_line_start"));
    entry.bits = static_cast<unsigned char>(bits[i]);
    entry.compressed_dictionary_size = (compressed_dictionary_size.size() == 0) ?
      0 :
      static_cast<std::uint16_t>(compressed_dictionary_size[i]);

    const size_t dictionary_offset = static_cast<size_t>(i) * fastqindex_core::WINDOW_SIZE;
    entry.dictionary.resize(fastqindex_core::WINDOW_SIZE);
    std::copy(
      dictionary_blob.begin() + static_cast<R_xlen_t>(dictionary_offset),
      dictionary_blob.begin() + static_cast<R_xlen_t>(dictionary_offset + fastqindex_core::WINDOW_SIZE),
      entry.dictionary.begin()
    );
    indexed_file.index.entries.push_back(std::move(entry));
  }

  indexed_file.index.total_lines = asUInt64(Rcpp::as<double>(payload["total_lines"]), "total_lines");
  indexed_file.record_count = asUInt64(Rcpp::as<double>(payload["record_count"]), "record_count");
  return indexed_file;
}

IndexBundle deserializeIndexBundle(const Rcpp::List& payload) {
  const int schema_version = Rcpp::as<int>(payload["schema_version"]);
  if (schema_version != 1 && schema_version != 2) {
    throw std::runtime_error("Unsupported serialized index schema version.");
  }

  IndexBundle bundle;
  bundle.format = Rcpp::as<std::string>(payload["format"]);
  bundle.record_size = Rcpp::as<int>(payload["record_size"]);
  bundle.n_records = asUInt64(Rcpp::as<double>(payload["n_records"]), "n_records");
  bundle.files = Rcpp::as<std::vector<std::string>>(payload["files"]);

  Rcpp::NumericVector record_offsets = payload["record_offsets"];
  bundle.record_offsets.reserve(static_cast<size_t>(record_offsets.size()));
  for (R_xlen_t i = 0; i < record_offsets.size(); ++i) {
    bundle.record_offsets.push_back(asUInt64(record_offsets[i], "record_offsets"));
  }

  Rcpp::List indexed_files = payload["indexed_files"];
  if (indexed_files.size() != static_cast<R_xlen_t>(bundle.files.size())) {
    throw std::runtime_error("Malformed payload: file and index counts differ.");
  }
  bundle.indexed_files.reserve(static_cast<size_t>(indexed_files.size()));
  for (R_xlen_t i = 0; i < indexed_files.size(); ++i) {
    bundle.indexed_files.push_back(deserializeIndexedFile(indexed_files[i], schema_version));
  }

  return bundle;
}

bool readGzLine(gzFile stream, std::string& out) {
  out.clear();
  char buffer[kGzBufferSize];
  while (true) {
    char* read = gzgets(stream, buffer, kGzBufferSize);
    if (read == nullptr) {
      return !out.empty();
    }
    out.append(buffer);
    if (!out.empty() && out.back() == '\n') {
      out.pop_back();
      if (!out.empty() && out.back() == '\r') {
        out.pop_back();
      }
      return true;
    }
  }
}

std::string trimPrefix(const std::string& line, char prefix) {
  if (!line.empty() && line.front() == prefix) {
    return line.substr(1);
  }
  return line;
}

std::string openAndDetectType(const std::string& path) {
  gzFile stream = gzopen(path.c_str(), "rb");
  if (stream == nullptr) {
    throw std::runtime_error("Could not open gzipped input: " + path);
  }
  std::string line;
  while (readGzLine(stream, line)) {
    if (line.empty()) {
      continue;
    }
    gzclose(stream);
    if (line.front() == '>') {
      return "fasta";
    }
    if (line.front() == '@') {
      return "fastq";
    }
    throw std::runtime_error("Cannot infer type from first non-empty line in: " + path);
  }
  gzclose(stream);
  throw std::runtime_error("Input appears empty: " + path);
}

Rcpp::DataFrame fileInfoDataFrame(const Rcpp::CharacterVector& files) {
  R_xlen_t n = files.size();
  Rcpp::CharacterVector path(n);
  Rcpp::NumericVector size(n);
  Rcpp::NumericVector mtime(n);
  Rcpp::Function file_size("file.size");
  Rcpp::Function file_mtime("file.mtime");
  for (R_xlen_t i = 0; i < n; ++i) {
    std::string p = Rcpp::as<std::string>(files[i]);
    path[i] = p;
    size[i] = Rcpp::as<double>(file_size(p));
    mtime[i] = Rcpp::as<double>(file_mtime(p));
  }
  return Rcpp::DataFrame::create(
    Rcpp::Named("path") = path,
    Rcpp::Named("size") = size,
    Rcpp::Named("mtime") = mtime
  );
}

IndexedFile readSingleFqi(const std::string& fqi_file) {
  auto source = fastqindex_core::FileSource::from(std::filesystem::path(fqi_file));
  fastqindex_core::IndexReader reader(source);
  std::vector<std::shared_ptr<fastqindex_core::IndexEntry>> upstream_entries =
    reader.readIndexFile();
  if (upstream_entries.empty()) {
    const auto errors = reader.getErrorMessages();
    std::string suffix = errors.empty() ? "" : (" (" + errors.front() + ")");
    throw std::runtime_error("Failed to parse .fqi index: " + fqi_file + suffix);
  }

  IndexedFile indexed_file;
  indexed_file.index.entries.reserve(upstream_entries.size());
  for (const auto& upstream_entry : upstream_entries) {
    indexed_file.index.entries.push_back(*upstream_entry);
  }
  indexed_file.index.total_lines = static_cast<std::uint64_t>(
    reader.getIndexHeader().linesInIndexedFile
  );
  return indexed_file;
}

std::uint64_t countGzLines(const std::string& path) {
  gzFile stream = gzopen(path.c_str(), "rb");
  if (stream == nullptr) {
    throw std::runtime_error("Could not open gzipped input: " + path);
  }
  std::uint64_t total_lines = 0;
  std::string line;
  while (readGzLine(stream, line)) {
    total_lines++;
  }
  gzclose(stream);
  return total_lines;
}

std::vector<std::pair<std::uint64_t, std::uint64_t>> buildDenseRegionsFromSortedUnique(
  const std::vector<std::uint64_t>& ids
) {
  std::vector<std::pair<std::uint64_t, std::uint64_t>> regions;
  if (ids.empty()) {
    return regions;
  }
  std::uint64_t start = ids[0];
  std::uint64_t last = ids[0];
  for (size_t i = 1; i < ids.size(); ++i) {
    if (ids[i] <= last + kDenseGap + 1) {
      last = ids[i];
    } else {
      regions.emplace_back(start, last);
      start = ids[i];
      last = ids[i];
    }
  }
  regions.emplace_back(start, last);
  return regions;
}

std::vector<std::pair<std::uint64_t, std::uint64_t>> denseRegions(std::vector<std::uint64_t> ids) {
  if (ids.empty()) {
    return {};
  }
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
  return buildDenseRegionsFromSortedUnique(ids);
}

ParsedRecord parseRecordFromFixedLines(
  const std::string& type,
  const std::vector<std::string>& lines,
  size_t offset,
  bool include_qual
) {
  if (type == "fastq") {
    if (offset + 3 >= lines.size()) {
      throw std::runtime_error("Malformed FASTQ lines in extracted region.");
    }
    if (include_qual) {
      return ParsedRecord{trimPrefix(lines[offset], '@'), lines[offset + 1], lines[offset + 3]};
    }
    return ParsedRecord{trimPrefix(lines[offset], '@'), lines[offset + 1], ""};
  }
  // fastqindexr change: FASTA extraction currently assumes one sequence line per record.
  if (offset + 1 >= lines.size()) {
    throw std::runtime_error("Malformed FASTA lines in extracted region.");
  }
  return ParsedRecord{trimPrefix(lines[offset], '>'), lines[offset + 1], ""};
}

void sequentialFallbackCollect(
  const std::string& file,
  const std::string& type,
  const std::unordered_set<std::uint64_t>& requested_local_ids,
  std::uint64_t global_offset,
  std::unordered_map<std::uint64_t, ParsedRecord>* selected_global,
  bool include_qual
) {
  gzFile stream = gzopen(file.c_str(), "rb");
  if (stream == nullptr) {
    throw std::runtime_error("Fallback parser could not open file: " + file);
  }

  std::uint64_t local_id = 0;
  std::string a;
  std::string b;
  std::string c;
  std::string d;

  if (type == "fastq") {
    while (readGzLine(stream, a)) {
      if (!readGzLine(stream, b) || !readGzLine(stream, c) || !readGzLine(stream, d)) {
        break;
      }
      if (requested_local_ids.find(local_id) != requested_local_ids.end()) {
        const std::string& q = include_qual ? d : "";
        (*selected_global)[global_offset + local_id] = ParsedRecord{trimPrefix(a, '@'), b, q};
      }
      local_id++;
    }
  } else {
    while (readGzLine(stream, a)) {
      if (!readGzLine(stream, b)) {
        break;
      }
      if (requested_local_ids.find(local_id) != requested_local_ids.end()) {
        (*selected_global)[global_offset + local_id] = ParsedRecord{trimPrefix(a, '>'), b, ""};
      }
      local_id++;
    }
  }

  gzclose(stream);
}

std::vector<std::uint64_t> parseRequestedIds(const Rcpp::NumericVector& ids_zero_based) {
  std::vector<std::uint64_t> requested(ids_zero_based.size());
  for (R_xlen_t i = 0; i < ids_zero_based.size(); ++i) {
    requested[static_cast<size_t>(i)] = static_cast<std::uint64_t>(ids_zero_based[i]);
  }
  return requested;
}

bool requestIsSortedUnique(const std::vector<std::uint64_t>& requested) {
  if (requested.size() <= 1) {
    return true;
  }
  if (!std::is_sorted(requested.begin(), requested.end())) {
    return false;
  }
  return std::adjacent_find(requested.begin(), requested.end()) == requested.end();
}

SelectedRecordMap collectRequestedRecords(
  const Rcpp::CharacterVector& files,
  const IndexBundle& bundle,
  const std::vector<std::uint64_t>& requested,
  bool include_qual
) {
  const bool need_qual = (bundle.format == "fastq" && include_qual);
  const bool sorted_unique = requestIsSortedUnique(requested);
  SelectedRecordMap selected_global;
  selected_global.reserve(requested.size());

  std::vector<std::vector<std::uint64_t>> local_ids_by_file(bundle.files.size());
  for (auto gid : requested) {
    auto up = std::upper_bound(bundle.record_offsets.begin(), bundle.record_offsets.end(), gid);
    if (up == bundle.record_offsets.begin() || up == bundle.record_offsets.end()) {
      throw std::runtime_error("Requested id out of range for indexed files.");
    }
    size_t file_idx = static_cast<size_t>((up - bundle.record_offsets.begin()) - 1);
    if (file_idx >= bundle.files.size()) {
      file_idx = bundle.files.size() - 1;
    }
    std::uint64_t local_id = gid - bundle.record_offsets[file_idx];
    local_ids_by_file[file_idx].push_back(local_id);
  }

  fastqindex_core::Extractor extractor;
  for (size_t file_idx = 0; file_idx < local_ids_by_file.size(); ++file_idx) {
    std::vector<std::uint64_t> local_ids = local_ids_by_file[file_idx];
    if (local_ids.empty()) {
      continue;
    }
    if (!sorted_unique) {
      std::sort(local_ids.begin(), local_ids.end());
      local_ids.erase(std::unique(local_ids.begin(), local_ids.end()), local_ids.end());
    }
    std::unordered_set<std::uint64_t> requested_set(local_ids.begin(), local_ids.end());
    auto regions = denseRegions(local_ids);
    for (const auto& rg : regions) {
      std::uint64_t region_start = rg.first;
      std::uint64_t region_end = rg.second;
      std::uint64_t line_start = region_start * static_cast<std::uint64_t>(bundle.record_size);
      std::uint64_t line_count =
        (region_end - region_start + 1) * static_cast<std::uint64_t>(bundle.record_size);

      try {
        std::vector<std::string> lines = extractor.extract(
          Rcpp::as<std::string>(files[file_idx]),
          bundle.indexed_files[file_idx].index.entries,
          line_start,
          line_count
        );
        if (lines.size() % static_cast<std::size_t>(bundle.record_size) != 0) {
          throw std::runtime_error("Malformed extracted line count for region.");
        }
        const std::size_t expected_records = static_cast<std::size_t>(
          region_end - region_start + 1
        );
        std::size_t extracted_records = lines.size() / static_cast<std::size_t>(bundle.record_size);
        if (extracted_records < expected_records) {
          throw std::runtime_error("Partial extraction for requested region.");
        }
        for (std::size_t r = 0; r < extracted_records; ++r) {
          std::uint64_t local_id = region_start + static_cast<std::uint64_t>(r);
          if (requested_set.find(local_id) == requested_set.end()) {
            continue;
          }
          ParsedRecord rec = parseRecordFromFixedLines(
            bundle.format,
            lines,
            r * static_cast<std::size_t>(bundle.record_size),
            need_qual
          );
          std::uint64_t global_id = bundle.record_offsets[file_idx] + local_id;
          selected_global[global_id] = std::move(rec);
        }
      } catch (...) {
        // fastqindexr change: Fallback preserves correctness when random-seek
        // extraction fails for particular gzip layouts.
        sequentialFallbackCollect(
          Rcpp::as<std::string>(files[file_idx]),
          bundle.format,
          requested_set,
          bundle.record_offsets[file_idx],
          &selected_global,
          need_qual
        );
      }
    }
  }

  return selected_global;
}

std::string resolveOutputType(const std::string& requested_type, const std::string& source_type) {
  std::string out_type = requested_type;
  if (out_type == "auto") {
    out_type = source_type;
  }
  if (out_type != "fasta" && out_type != "fastq") {
    throw std::runtime_error("`type` must be one of 'auto', 'fasta', or 'fastq'.");
  }
  if (source_type == "fasta" && out_type == "fastq") {
    throw std::runtime_error("Cannot emit FASTQ output from FASTA index input.");
  }
  return out_type;
}

std::string renderRecord(const ParsedRecord& record, const std::string& output_type) {
  if (output_type == "fastq") {
    return "@" + record.id + "\n" + record.seq + "\n+\n" + record.qual + "\n";
  }
  return ">" + record.id + "\n" + record.seq + "\n";
}

std::unique_ptr<OutputWriter> makeOutputWriter(
  const std::string& path,
  bool append,
  bool compress
) {
  if (compress) {
    return std::make_unique<GzipOutputWriter>(path, append);
  }
  return std::make_unique<PlainOutputWriter>(path, append);
}

std::uint64_t sequentialFallbackStreamToFile(
  const std::string& file,
  const std::string& type,
  const std::unordered_set<std::uint64_t>& requested_local_ids,
  bool need_qual,
  const std::string& output_type,
  OutputWriter* writer
) {
  std::uint64_t count = 0;
  gzFile stream = gzopen(file.c_str(), "rb");
  if (stream == nullptr) {
    throw std::runtime_error("Fallback stream could not open file: " + file);
  }

  std::uint64_t local_id = 0;
  std::string a;
  std::string b;
  std::string c;
  std::string d;

  if (type == "fastq") {
    while (readGzLine(stream, a)) {
      if (!readGzLine(stream, b) || !readGzLine(stream, c) || !readGzLine(stream, d)) {
        break;
      }
      if (requested_local_ids.find(local_id) != requested_local_ids.end()) {
        const std::string& q = need_qual ? d : "";
        const ParsedRecord rec{trimPrefix(a, '@'), b, q};
        writer->write(renderRecord(rec, output_type));
        count++;
      }
      local_id++;
    }
  } else {
    while (readGzLine(stream, a)) {
      if (!readGzLine(stream, b)) {
        break;
      }
      if (requested_local_ids.find(local_id) != requested_local_ids.end()) {
        const ParsedRecord rec{trimPrefix(a, '>'), b, ""};
        writer->write(renderRecord(rec, output_type));
        count++;
      }
      local_id++;
    }
  }

  gzclose(stream);
  return count;
}

double fullCollectAndRewriteToFile(
  Rcpp::CharacterVector files,
  const IndexBundle& bundle,
  const std::vector<std::uint64_t>& requested,
  const std::string& resolved_output_type,
  const std::string& outfile,
  bool compress
) {
  const bool include_qual = (resolved_output_type == "fastq");
  const SelectedRecordMap selected_global = collectRequestedRecords(
    files, bundle, requested, include_qual
  );
  std::unique_ptr<OutputWriter> writer = makeOutputWriter(outfile, false, compress);
  std::uint64_t written = 0;
  for (auto gid : requested) {
    auto it = selected_global.find(gid);
    if (it == selected_global.end()) {
      throw std::runtime_error("Could not resolve all requested ids with provided index/files.");
    }
    writer->write(renderRecord(it->second, resolved_output_type));
    written++;
  }
  return static_cast<double>(written);
}

double streamSortedUniqueToFile(
  Rcpp::CharacterVector files,
  const IndexBundle& bundle,
  const std::vector<std::uint64_t>& requested,
  const std::string& resolved_output_type,
  const std::string& outfile,
  bool compress
) {
  const bool include_qual = (resolved_output_type == "fastq");
  const bool need_qual = (bundle.format == "fastq" && include_qual);
  std::uint64_t n_written = 0;
  std::unique_ptr<OutputWriter> writer = makeOutputWriter(outfile, false, compress);

  std::vector<std::vector<std::uint64_t>> local_ids_by_file(bundle.files.size());
  for (auto gid : requested) {
    auto up = std::upper_bound(bundle.record_offsets.begin(), bundle.record_offsets.end(), gid);
    if (up == bundle.record_offsets.begin() || up == bundle.record_offsets.end()) {
      throw std::runtime_error("Requested id out of range for indexed files.");
    }
    size_t file_idx = static_cast<size_t>((up - bundle.record_offsets.begin()) - 1);
    if (file_idx >= bundle.files.size()) {
      file_idx = bundle.files.size() - 1;
    }
    std::uint64_t local_id = gid - bundle.record_offsets[file_idx];
    local_ids_by_file[file_idx].push_back(local_id);
  }

  fastqindex_core::Extractor extractor;
  for (size_t file_idx = 0; file_idx < local_ids_by_file.size(); ++file_idx) {
    std::vector<std::uint64_t> local_ids = std::move(local_ids_by_file[file_idx]);
    if (local_ids.empty()) {
      continue;
    }
    const std::unordered_set<std::uint64_t> requested_set(local_ids.begin(), local_ids.end());
    const auto regions = buildDenseRegionsFromSortedUnique(local_ids);
    bool wrote_in_this_file = false;

    for (const auto& rg : regions) {
      const std::uint64_t region_start = rg.first;
      const std::uint64_t region_end = rg.second;
      const std::uint64_t line_start = region_start * static_cast<std::uint64_t>(bundle.record_size);
      const std::uint64_t line_count =
        (region_end - region_start + 1) * static_cast<std::uint64_t>(bundle.record_size);
      try {
        std::vector<std::string> lines = extractor.extract(
          Rcpp::as<std::string>(files[file_idx]),
          bundle.indexed_files[file_idx].index.entries,
          line_start,
          line_count
        );
        if (lines.size() % static_cast<std::size_t>(bundle.record_size) != 0) {
          throw std::runtime_error("Malformed extracted line count for region.");
        }
        const std::size_t expected_records = static_cast<std::size_t>(
          region_end - region_start + 1
        );
        const std::size_t extracted_records =
          lines.size() / static_cast<std::size_t>(bundle.record_size);
        if (extracted_records < expected_records) {
          throw std::runtime_error("Partial extraction for requested region.");
        }
        for (std::size_t r = 0; r < extracted_records; ++r) {
          const std::uint64_t rec_local = region_start + static_cast<std::uint64_t>(r);
          if (requested_set.find(rec_local) == requested_set.end()) {
            continue;
          }
          const ParsedRecord rec = parseRecordFromFixedLines(
            bundle.format,
            lines,
            r * static_cast<std::size_t>(bundle.record_size),
            need_qual
          );
          writer->write(renderRecord(rec, resolved_output_type));
          wrote_in_this_file = true;
          n_written++;
        }
      } catch (...) {
        if (!wrote_in_this_file) {
          n_written += sequentialFallbackStreamToFile(
            Rcpp::as<std::string>(files[file_idx]),
            bundle.format,
            requested_set,
            need_qual,
            resolved_output_type,
            writer.get()
          );
        } else {
          return fullCollectAndRewriteToFile(
            files, bundle, requested, resolved_output_type, outfile, compress
          );
        }
        break;
      }
    }
  }

  return static_cast<double>(n_written);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_create_index(Rcpp::CharacterVector files, std::string type) {
  if (files.size() < 1) {
    throw std::runtime_error("`files` must contain at least one path.");
  }

  std::vector<std::string> paths(files.size());
  for (R_xlen_t i = 0; i < files.size(); ++i) {
    paths[i] = Rcpp::as<std::string>(files[i]);
  }

  if (type == "auto") {
    type = openAndDetectType(paths.front());
  }

  int record_size = (type == "fastq") ? 4 : 2;

  fastqindex_core::Indexer indexer;
  IndexBundle bundle;
  bundle.format = type;
  bundle.record_size = record_size;
  bundle.files = paths;
  bundle.record_offsets.push_back(0);

  Rcpp::NumericVector per_file_counts(paths.size());
  for (size_t i = 0; i < paths.size(); ++i) {
    IndexedFile idxf;
    idxf.index = indexer.createIndex(paths[i]);
    if (idxf.index.total_lines % static_cast<std::uint64_t>(record_size) != 0) {
      throw std::runtime_error("Indexed line count is not a multiple of record size for: " + paths[i]);
    }
    idxf.record_count = idxf.index.total_lines / static_cast<std::uint64_t>(record_size);
    bundle.indexed_files.push_back(std::move(idxf));
    bundle.n_records += bundle.indexed_files.back().record_count;
    bundle.record_offsets.push_back(bundle.n_records);
    per_file_counts[i] = static_cast<double>(bundle.indexed_files.back().record_count);
  }

  Rcpp::XPtr<IndexBundle> index_ptr(new IndexBundle(std::move(bundle)), true);
  Rcpp::List payload = serializeIndexBundle(*index_ptr);

  Rcpp::NumericVector offsets(per_file_counts.size() + 1);
  for (R_xlen_t i = 0; i < offsets.size(); ++i) {
    offsets[i] = static_cast<double>((*index_ptr).record_offsets[static_cast<size_t>(i)]);
  }

  return Rcpp::List::create(
    Rcpp::Named("format") = type,
    Rcpp::Named("files") = files,
    Rcpp::Named("file_info") = fileInfoDataFrame(files),
    Rcpp::Named("file_record_counts") = per_file_counts,
    Rcpp::Named("file_record_offsets") = offsets,
    Rcpp::Named("n_records") = static_cast<double>((*index_ptr).n_records),
    Rcpp::Named("index_payload") = payload,
    Rcpp::Named("index_ptr") = index_ptr,
    Rcpp::Named("record_size") = record_size,
    Rcpp::Named("engine") = "upstream_adapted_v2"
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_read_fqi_index(
  Rcpp::CharacterVector fqi_files,
  Rcpp::CharacterVector files,
  std::string type
) {
  if (fqi_files.size() < 1) {
    throw std::runtime_error("`fqi_files` must contain at least one path.");
  }
  if (files.size() != fqi_files.size()) {
    throw std::runtime_error("`files` must have the same length as `fqi_files`.");
  }

  std::vector<std::string> fqi_paths(fqi_files.size());
  std::vector<std::string> data_paths(files.size());
  for (R_xlen_t i = 0; i < fqi_files.size(); ++i) {
    fqi_paths[i] = Rcpp::as<std::string>(fqi_files[i]);
    data_paths[i] = Rcpp::as<std::string>(files[i]);
  }

  if (type != "fasta" && type != "fastq") {
    throw std::runtime_error("`type` must resolve to 'fasta' or 'fastq'.");
  }
  const int record_size = (type == "fastq") ? 4 : 2;

  IndexBundle bundle;
  bundle.format = type;
  bundle.record_size = record_size;
  bundle.files = data_paths;
  bundle.record_offsets.push_back(0);

  Rcpp::NumericVector per_file_counts(data_paths.size());
  for (size_t i = 0; i < fqi_paths.size(); ++i) {
    IndexedFile idxf = readSingleFqi(fqi_paths[i]);
    if (idxf.index.total_lines == 0 && std::filesystem::exists(data_paths[i])) {
      idxf.index.total_lines = countGzLines(data_paths[i]);
    }
    if (idxf.index.total_lines > 0 &&
        idxf.index.total_lines % static_cast<std::uint64_t>(record_size) != 0) {
      throw std::runtime_error(
        "Indexed line count is not a multiple of record size for: " + fqi_paths[i]
      );
    }
    idxf.record_count = idxf.index.total_lines / static_cast<std::uint64_t>(record_size);
    bundle.indexed_files.push_back(std::move(idxf));
    bundle.n_records += bundle.indexed_files.back().record_count;
    bundle.record_offsets.push_back(bundle.n_records);
    per_file_counts[i] = static_cast<double>(bundle.indexed_files.back().record_count);
  }

  Rcpp::XPtr<IndexBundle> index_ptr(new IndexBundle(std::move(bundle)), true);
  Rcpp::List payload = serializeIndexBundle(*index_ptr);

  Rcpp::NumericVector offsets(per_file_counts.size() + 1);
  for (R_xlen_t i = 0; i < offsets.size(); ++i) {
    offsets[i] = static_cast<double>((*index_ptr).record_offsets[static_cast<size_t>(i)]);
  }

  return Rcpp::List::create(
    Rcpp::Named("format") = type,
    Rcpp::Named("files") = files,
    Rcpp::Named("fqi_files") = fqi_files,
    Rcpp::Named("file_info") = fileInfoDataFrame(files),
    Rcpp::Named("file_record_counts") = per_file_counts,
    Rcpp::Named("file_record_offsets") = offsets,
    Rcpp::Named("n_records") = static_cast<double>((*index_ptr).n_records),
    Rcpp::Named("index_payload") = payload,
    Rcpp::Named("index_ptr") = index_ptr,
    Rcpp::Named("record_size") = record_size,
    Rcpp::Named("engine") = "upstream_cli_fqi_v1"
  );
}

// [[Rcpp::export]]
SEXP cpp_restore_index_ptr(Rcpp::List index_payload) {
  IndexBundle bundle = deserializeIndexBundle(index_payload);
  Rcpp::XPtr<IndexBundle> index_ptr(new IndexBundle(std::move(bundle)), true);
  return index_ptr;
}

// [[Rcpp::export]]
bool cpp_index_ptr_is_valid(SEXP index_ptr) {
  return TYPEOF(index_ptr) == EXTPTRSXP && R_ExternalPtrAddr(index_ptr) != nullptr;
}

// [[Rcpp::export]]
Rcpp::List cpp_extract_sequences(
  Rcpp::CharacterVector files,
  std::string type,
  Rcpp::NumericVector ids_zero_based,
  SEXP index_ptr_sexp,
  bool include_qual
) {
  if (!cpp_index_ptr_is_valid(index_ptr_sexp)) {
    throw std::runtime_error("Invalid index pointer; recreate from serialized payload.");
  }
  Rcpp::XPtr<IndexBundle> index_ptr(index_ptr_sexp);
  IndexBundle& bundle = *index_ptr;
  if (bundle.format != type) {
    throw std::runtime_error("Type mismatch between index and request.");
  }
  if (static_cast<size_t>(files.size()) != bundle.files.size()) {
    throw std::runtime_error("Override `file` must contain same file count as index.");
  }

  const bool return_qual = (bundle.format == "fastq" && include_qual);
  if (ids_zero_based.size() < 1) {
    if (return_qual) {
      return Rcpp::List::create(
        Rcpp::Named("seq_id") = Rcpp::CharacterVector(),
        Rcpp::Named("seq") = Rcpp::CharacterVector(),
        Rcpp::Named("qual") = Rcpp::CharacterVector()
      );
    }
    return Rcpp::List::create(
      Rcpp::Named("seq_id") = Rcpp::CharacterVector(),
      Rcpp::Named("seq") = Rcpp::CharacterVector()
    );
  }

  const std::vector<std::uint64_t> requested = parseRequestedIds(ids_zero_based);
  const bool need_parsed_qual = (bundle.format == "fastq" && include_qual);
  const SelectedRecordMap selected_global = collectRequestedRecords(
    files, bundle, requested, need_parsed_qual
  );

  R_xlen_t n = ids_zero_based.size();
  Rcpp::CharacterVector out_id(n);
  Rcpp::CharacterVector out_seq(n);
  Rcpp::CharacterVector out_qual;
  if (return_qual) {
    out_qual = Rcpp::CharacterVector(n);
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    std::uint64_t gid = requested[i];
    auto it = selected_global.find(gid);
    if (it == selected_global.end()) {
      throw std::runtime_error("Could not resolve all requested ids with provided index/files.");
    }
    out_id[i] = it->second.id;
    out_seq[i] = it->second.seq;
    if (return_qual) {
      out_qual[i] = it->second.qual;
    }
  }

  if (return_qual) {
    return Rcpp::List::create(
      Rcpp::Named("seq_id") = out_id,
      Rcpp::Named("seq") = out_seq,
      Rcpp::Named("qual") = out_qual
    );
  }
  return Rcpp::List::create(
    Rcpp::Named("seq_id") = out_id,
    Rcpp::Named("seq") = out_seq
  );
}

// [[Rcpp::export]]
double cpp_extract_sequences_to_file(
  Rcpp::CharacterVector files,
  std::string source_type,
  Rcpp::NumericVector ids_zero_based,
  SEXP index_ptr_sexp,
  std::string output_type,
  std::string outfile,
  bool append,
  bool compress
) {
  if (!cpp_index_ptr_is_valid(index_ptr_sexp)) {
    throw std::runtime_error("Invalid index pointer; recreate from serialized payload.");
  }
  Rcpp::XPtr<IndexBundle> index_ptr(index_ptr_sexp);
  IndexBundle& bundle = *index_ptr;
  if (bundle.format != source_type) {
    throw std::runtime_error("Type mismatch between index and request.");
  }
  if (static_cast<size_t>(files.size()) != bundle.files.size()) {
    throw std::runtime_error("Override `file` must contain same file count as index.");
  }
  if (ids_zero_based.size() < 1) {
    return 0.0;
  }

  const std::string resolved_output_type = resolveOutputType(output_type, bundle.format);
  const std::vector<std::uint64_t> requested = parseRequestedIds(ids_zero_based);
  const bool need_parsed_qual = (resolved_output_type == "fastq");

  if (requestIsSortedUnique(requested) && !append) {
    return streamSortedUniqueToFile(
      files, bundle, requested, resolved_output_type, outfile, compress
    );
  }

  const SelectedRecordMap selected_global = collectRequestedRecords(
    files, bundle, requested, need_parsed_qual
  );

  std::unique_ptr<OutputWriter> writer;
  if (compress) {
    writer = std::make_unique<GzipOutputWriter>(outfile, append);
  } else {
    writer = std::make_unique<PlainOutputWriter>(outfile, append);
  }

  std::uint64_t written = 0;
  for (auto gid : requested) {
    auto it = selected_global.find(gid);
    if (it == selected_global.end()) {
      throw std::runtime_error("Could not resolve all requested ids with provided index/files.");
    }
    writer->write(renderRecord(it->second, resolved_output_type));
    written++;
  }
  return static_cast<double>(written);
}
