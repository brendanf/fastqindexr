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

enum class ExtractExecutionMode {
  Indexed = 0,
  SequentialOnly = 1
};

struct ExtractionDiagnostics {
  std::uint64_t regions_planned{0};
  std::uint64_t extract_attempts{0};
  std::uint64_t extract_failures{0};
  std::uint64_t extract_failures_partial{0};
  std::uint64_t extract_failures_malformed{0};
  std::uint64_t extract_failures_other{0};
  std::uint64_t fallback_invocations{0};
  std::uint64_t fallback_records{0};
  std::string last_failure_message{};
};

void recordExtractFailure(ExtractionDiagnostics* diag, const std::string& msg) {
  if (diag == nullptr) {
    return;
  }
  diag->extract_failures++;
  diag->last_failure_message = msg;
  if (msg.find("Partial extraction for requested region.") != std::string::npos) {
    diag->extract_failures_partial++;
  } else if (msg.find("Malformed extracted line count for region.") != std::string::npos) {
    diag->extract_failures_malformed++;
  } else {
    diag->extract_failures_other++;
  }
}

Rcpp::List diagnosticsToList(const ExtractionDiagnostics& diag) {
  return Rcpp::List::create(
    Rcpp::Named("regions_planned") = static_cast<double>(diag.regions_planned),
    Rcpp::Named("extract_attempts") = static_cast<double>(diag.extract_attempts),
    Rcpp::Named("extract_failures") = static_cast<double>(diag.extract_failures),
    Rcpp::Named("extract_failures_partial") = static_cast<double>(diag.extract_failures_partial),
    Rcpp::Named("extract_failures_malformed") = static_cast<double>(diag.extract_failures_malformed),
    Rcpp::Named("extract_failures_other") = static_cast<double>(diag.extract_failures_other),
    Rcpp::Named("fallback_invocations") = static_cast<double>(diag.fallback_invocations),
    Rcpp::Named("fallback_records") = static_cast<double>(diag.fallback_records),
    Rcpp::Named("last_failure_message") = diag.last_failure_message
  );
}

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
  const std::vector<std::uint64_t>& ids,
  std::uint64_t max_bridge_gap,
  std::uint64_t max_region_bytes,
  int record_size,
  ExtractionDiagnostics* diag
) {
  std::vector<std::pair<std::uint64_t, std::uint64_t>> regions;
  if (ids.empty()) {
    return regions;
  }
  if (record_size <= 0) {
    throw std::runtime_error("Invalid record size for region planning.");
  }
  std::uint64_t start = ids[0];
  std::uint64_t last = ids[0];
  const std::uint64_t bytes_per_record = static_cast<std::uint64_t>(record_size);
  for (size_t i = 1; i < ids.size(); ++i) {
    const bool within_gap = (ids[i] <= last + max_bridge_gap + 1);
    const std::uint64_t span_records = ids[i] - start + 1;
    const bool within_bytes = (span_records <= max_region_bytes / bytes_per_record);
    if (within_gap && within_bytes) {
      last = ids[i];
    } else {
      regions.emplace_back(start, last);
      start = ids[i];
      last = ids[i];
    }
  }
  regions.emplace_back(start, last);
  if (diag != nullptr) {
    diag->regions_planned += static_cast<std::uint64_t>(regions.size());
  }
  return regions;
}

std::vector<std::pair<std::uint64_t, std::uint64_t>> denseRegions(
  std::vector<std::uint64_t> ids,
  std::uint64_t max_bridge_gap,
  std::uint64_t max_region_bytes,
  int record_size,
  ExtractionDiagnostics* diag
) {
  if (ids.empty()) {
    return {};
  }
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
  return buildDenseRegionsFromSortedUnique(
    ids,
    max_bridge_gap,
    max_region_bytes,
    record_size,
    diag
  );
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
  bool include_qual,
  ExtractionDiagnostics* diag
) {
  if (diag != nullptr) {
    diag->fallback_invocations++;
  }
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
        if (diag != nullptr) {
          diag->fallback_records++;
        }
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
        if (diag != nullptr) {
          diag->fallback_records++;
        }
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

std::uint64_t parseNonNegativeWhole(double value, const std::string& field) {
  if (!std::isfinite(value) || value < 0.0 || std::floor(value) != value) {
    throw std::runtime_error("`" + field + "` must be a non-negative whole number.");
  }
  return static_cast<std::uint64_t>(value);
}

std::uint64_t parsePositiveWhole(double value, const std::string& field) {
  if (!std::isfinite(value) || value <= 0.0 || std::floor(value) != value) {
    throw std::runtime_error("`" + field + "` must be a positive whole number.");
  }
  return static_cast<std::uint64_t>(value);
}

ExtractExecutionMode parseExtractExecutionMode(const std::string& mode) {
  if (mode == "indexed") {
    return ExtractExecutionMode::Indexed;
  }
  if (mode == "sequential_only") {
    return ExtractExecutionMode::SequentialOnly;
  }
  throw std::runtime_error("`extract_mode` must be 'indexed' or 'sequential_only'.");
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
  bool include_qual,
  std::uint64_t max_bridge_gap,
  std::uint64_t max_region_bytes,
  ExtractExecutionMode mode,
  ExtractionDiagnostics* diag
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
    if (mode == ExtractExecutionMode::SequentialOnly) {
      sequentialFallbackCollect(
        Rcpp::as<std::string>(files[file_idx]),
        bundle.format,
        requested_set,
        bundle.record_offsets[file_idx],
        &selected_global,
        need_qual,
        diag
      );
      continue;
    }
    auto regions = denseRegions(
      local_ids,
      max_bridge_gap,
      max_region_bytes,
      bundle.record_size,
      diag
    );
    for (const auto& rg : regions) {
      std::uint64_t region_start = rg.first;
      std::uint64_t region_end = rg.second;
      std::uint64_t line_start = region_start * static_cast<std::uint64_t>(bundle.record_size);
      std::uint64_t line_count =
        (region_end - region_start + 1) * static_cast<std::uint64_t>(bundle.record_size);

      try {
        if (diag != nullptr) {
          diag->extract_attempts++;
        }
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
      } catch (const std::exception& e) {
        recordExtractFailure(diag, e.what());
        // fastqindexr change: Fallback preserves correctness when random-seek
        // extraction fails for particular gzip layouts.
        sequentialFallbackCollect(
          Rcpp::as<std::string>(files[file_idx]),
          bundle.format,
          requested_set,
          bundle.record_offsets[file_idx],
          &selected_global,
          need_qual,
          diag
        );
        break;
      } catch (...) {
        recordExtractFailure(diag, "Unknown indexed extraction failure.");
        sequentialFallbackCollect(
          Rcpp::as<std::string>(files[file_idx]),
          bundle.format,
          requested_set,
          bundle.record_offsets[file_idx],
          &selected_global,
          need_qual,
          diag
        );
        break;
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
  OutputWriter* writer,
  ExtractionDiagnostics* diag
) {
  if (diag != nullptr) {
    diag->fallback_invocations++;
  }
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
        if (diag != nullptr) {
          diag->fallback_records++;
        }
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
        if (diag != nullptr) {
          diag->fallback_records++;
        }
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
  bool compress,
  std::uint64_t max_bridge_gap,
  std::uint64_t max_region_bytes,
  ExtractExecutionMode mode,
  ExtractionDiagnostics* diag
) {
  const bool include_qual = (resolved_output_type == "fastq");
  const SelectedRecordMap selected_global = collectRequestedRecords(
    files,
    bundle,
    requested,
    include_qual,
    max_bridge_gap,
    max_region_bytes,
    mode,
    diag
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
  bool compress,
  std::uint64_t max_bridge_gap,
  std::uint64_t max_region_bytes,
  ExtractExecutionMode mode,
  ExtractionDiagnostics* diag
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
    if (mode == ExtractExecutionMode::SequentialOnly) {
      n_written += sequentialFallbackStreamToFile(
        Rcpp::as<std::string>(files[file_idx]),
        bundle.format,
        requested_set,
        need_qual,
        resolved_output_type,
        writer.get(),
        diag
      );
      continue;
    }
    const auto regions = buildDenseRegionsFromSortedUnique(
      local_ids,
      max_bridge_gap,
      max_region_bytes,
      bundle.record_size,
      diag
    );
    bool wrote_in_this_file = false;

    for (const auto& rg : regions) {
      const std::uint64_t region_start = rg.first;
      const std::uint64_t region_end = rg.second;
      const std::uint64_t line_start = region_start * static_cast<std::uint64_t>(bundle.record_size);
      const std::uint64_t line_count =
        (region_end - region_start + 1) * static_cast<std::uint64_t>(bundle.record_size);
      try {
        if (diag != nullptr) {
          diag->extract_attempts++;
        }
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
      } catch (const std::exception& e) {
        recordExtractFailure(diag, e.what());
        if (!wrote_in_this_file) {
          n_written += sequentialFallbackStreamToFile(
            Rcpp::as<std::string>(files[file_idx]),
            bundle.format,
            requested_set,
            need_qual,
            resolved_output_type,
            writer.get(),
            diag
          );
        } else {
          return fullCollectAndRewriteToFile(
            files,
            bundle,
            requested,
            resolved_output_type,
            outfile,
            compress,
            max_bridge_gap,
            max_region_bytes,
            mode,
            diag
          );
        }
        break;
      } catch (...) {
        recordExtractFailure(diag, "Unknown indexed extraction failure.");
        if (!wrote_in_this_file) {
          n_written += sequentialFallbackStreamToFile(
            Rcpp::as<std::string>(files[file_idx]),
            bundle.format,
            requested_set,
            need_qual,
            resolved_output_type,
            writer.get(),
            diag
          );
        } else {
          return fullCollectAndRewriteToFile(
            files,
            bundle,
            requested,
            resolved_output_type,
            outfile,
            compress,
            max_bridge_gap,
            max_region_bytes,
            mode,
            diag
          );
        }
        break;
      }
    }
  }

  return static_cast<double>(n_written);
}

SEXP makeDNAStringSetChunk(
  const std::vector<std::string>& seq_chunk,
  const std::vector<std::string>& id_chunk,
  const std::string& renumber_mode,
  std::uint64_t output_index_start
) {
  const R_xlen_t n = static_cast<R_xlen_t>(seq_chunk.size());
  Rcpp::CharacterVector seq_vec(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    seq_vec[i] = seq_chunk[static_cast<size_t>(i)];
  }

  Rcpp::Environment biostrings = Rcpp::Environment::namespace_env("Biostrings");
  Rcpp::Function dna_string_set = biostrings["DNAStringSet"];
  SEXP dna = dna_string_set(seq_vec, Rcpp::Named("use.names") = true);

  Rcpp::CharacterVector out_names(n);
  if (renumber_mode == "none") {
    for (R_xlen_t i = 0; i < n; ++i) {
      out_names[i] = id_chunk[static_cast<size_t>(i)];
    }
  } else if (renumber_mode == "zero_based") {
    for (R_xlen_t i = 0; i < n; ++i) {
      out_names[i] = std::to_string(output_index_start + static_cast<std::uint64_t>(i));
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      out_names[i] = std::to_string(output_index_start + static_cast<std::uint64_t>(i) + 1);
    }
  }
  Rcpp::Function set_names = Rcpp::Function("names<-");
  return set_names(dna, out_names);
}

void flushDNAChunk(
  std::vector<std::string>* seq_chunk,
  std::vector<std::string>* id_chunk,
  Rcpp::List* dna_chunks,
  const std::string& renumber_mode,
  std::uint64_t output_index_start
) {
  if (seq_chunk->empty()) {
    return;
  }
  dna_chunks->push_back(
    makeDNAStringSetChunk(*seq_chunk, *id_chunk, renumber_mode, output_index_start)
  );
  seq_chunk->clear();
  id_chunk->clear();
}

void appendRecordToDNAChunks(
  const ParsedRecord& rec,
  std::uint64_t chunk_chars_limit,
  std::vector<std::string>* seq_chunk,
  std::vector<std::string>* id_chunk,
  std::uint64_t* chunk_chars,
  std::uint64_t* output_index,
  Rcpp::List* dna_chunks,
  const std::string& renumber_mode
) {
  const std::uint64_t rec_chars = static_cast<std::uint64_t>(rec.seq.size());
  if (!seq_chunk->empty() && (*chunk_chars + rec_chars > chunk_chars_limit)) {
    const std::uint64_t chunk_start = *output_index - seq_chunk->size();
    flushDNAChunk(seq_chunk, id_chunk, dna_chunks, renumber_mode, chunk_start);
    *chunk_chars = 0;
  }

  seq_chunk->push_back(rec.seq);
  id_chunk->push_back(rec.id);
  *chunk_chars += rec_chars;
  (*output_index)++;
}

SEXP buildDNAStringSetFromRequested(
  Rcpp::CharacterVector files,
  const IndexBundle& bundle,
  const std::vector<std::uint64_t>& requested,
  std::uint64_t chunk_chars_limit,
  const std::string& renumber_mode,
  std::uint64_t max_bridge_gap,
  std::uint64_t max_region_bytes,
  ExtractExecutionMode mode,
  ExtractionDiagnostics* diag
) {
  const bool sorted_unique = requestIsSortedUnique(requested);
  std::vector<std::string> seq_chunk;
  std::vector<std::string> id_chunk;
  seq_chunk.reserve(1024);
  id_chunk.reserve(1024);
  Rcpp::List dna_chunks;
  std::uint64_t chunk_chars = 0;
  std::uint64_t output_index = 0;

  if (!sorted_unique) {
    const SelectedRecordMap selected_global = collectRequestedRecords(
      files,
      bundle,
      requested,
      false,
      max_bridge_gap,
      max_region_bytes,
      mode,
      diag
    );
    for (auto gid : requested) {
      auto it = selected_global.find(gid);
      if (it == selected_global.end()) {
        throw std::runtime_error("Could not resolve all requested ids with provided index/files.");
      }
      appendRecordToDNAChunks(
        it->second,
        chunk_chars_limit,
        &seq_chunk,
        &id_chunk,
        &chunk_chars,
        &output_index,
        &dna_chunks,
        renumber_mode
      );
    }
  } else {
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
      const std::uint64_t local_id = gid - bundle.record_offsets[file_idx];
      local_ids_by_file[file_idx].push_back(local_id);
    }

    fastqindex_core::Extractor extractor;
    for (size_t file_idx = 0; file_idx < local_ids_by_file.size(); ++file_idx) {
      const std::vector<std::uint64_t>& local_ids = local_ids_by_file[file_idx];
      if (local_ids.empty()) {
        continue;
      }
      const std::unordered_set<std::uint64_t> requested_set(local_ids.begin(), local_ids.end());
      if (mode == ExtractExecutionMode::SequentialOnly) {
        std::unordered_map<std::uint64_t, ParsedRecord> selected_global;
        sequentialFallbackCollect(
          Rcpp::as<std::string>(files[file_idx]),
          bundle.format,
          requested_set,
          bundle.record_offsets[file_idx],
          &selected_global,
          false,
          diag
        );
        for (auto rec_local : local_ids) {
          const std::uint64_t gid = bundle.record_offsets[file_idx] + rec_local;
          auto it = selected_global.find(gid);
          if (it == selected_global.end()) {
            throw std::runtime_error("Could not resolve all requested ids with provided index/files.");
          }
          appendRecordToDNAChunks(
            it->second,
            chunk_chars_limit,
            &seq_chunk,
            &id_chunk,
            &chunk_chars,
            &output_index,
            &dna_chunks,
            renumber_mode
          );
        }
        continue;
      }
      const auto regions = buildDenseRegionsFromSortedUnique(
        local_ids,
        max_bridge_gap,
        max_region_bytes,
        bundle.record_size,
        diag
      );
      bool fell_back = false;

      for (const auto& rg : regions) {
        const std::uint64_t region_start = rg.first;
        const std::uint64_t region_end = rg.second;
        const std::uint64_t line_start = region_start * static_cast<std::uint64_t>(bundle.record_size);
        const std::uint64_t line_count =
          (region_end - region_start + 1) * static_cast<std::uint64_t>(bundle.record_size);
        try {
          if (diag != nullptr) {
            diag->extract_attempts++;
          }
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
              false
            );
            appendRecordToDNAChunks(
              rec,
              chunk_chars_limit,
              &seq_chunk,
              &id_chunk,
              &chunk_chars,
              &output_index,
              &dna_chunks,
              renumber_mode
            );
          }
        } catch (const std::exception& e) {
          recordExtractFailure(diag, e.what());
          fell_back = true;
          break;
        } catch (...) {
          recordExtractFailure(diag, "Unknown indexed extraction failure.");
          fell_back = true;
          break;
        }
      }

      if (fell_back) {
        std::unordered_map<std::uint64_t, ParsedRecord> selected_global;
        sequentialFallbackCollect(
          Rcpp::as<std::string>(files[file_idx]),
          bundle.format,
          requested_set,
          bundle.record_offsets[file_idx],
          &selected_global,
          false,
          diag
        );
        for (auto rec_local : local_ids) {
          const std::uint64_t gid = bundle.record_offsets[file_idx] + rec_local;
          auto it = selected_global.find(gid);
          if (it == selected_global.end()) {
            throw std::runtime_error("Could not resolve all requested ids with provided index/files.");
          }
          appendRecordToDNAChunks(
            it->second,
            chunk_chars_limit,
            &seq_chunk,
            &id_chunk,
            &chunk_chars,
            &output_index,
            &dna_chunks,
            renumber_mode
          );
        }
      }
    }
  }

  if (!seq_chunk.empty()) {
    const std::uint64_t chunk_start = output_index - seq_chunk.size();
    flushDNAChunk(&seq_chunk, &id_chunk, &dna_chunks, renumber_mode, chunk_start);
  }

  if (dna_chunks.size() == 0) {
    Rcpp::Environment biostrings = Rcpp::Environment::namespace_env("Biostrings");
    Rcpp::Function dna_string_set = biostrings["DNAStringSet"];
    return dna_string_set(Rcpp::CharacterVector(), Rcpp::Named("use.names") = true);
  }
  if (dna_chunks.size() == 1) {
    return dna_chunks[0];
  }
  Rcpp::Environment base = Rcpp::Environment::base_env();
  Rcpp::Function do_call = base["do.call"];
  Rcpp::Function concat = base["c"];
  return do_call(concat, dna_chunks);
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
  bool include_qual,
  double max_bridge_gap,
  double max_region_bytes,
  std::string extract_mode,
  bool diagnostics
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
  const std::uint64_t bridge_gap = parseNonNegativeWhole(max_bridge_gap, "max_bridge_gap");
  const std::uint64_t region_bytes = parsePositiveWhole(max_region_bytes, "max_region_bytes");
  const ExtractExecutionMode mode = parseExtractExecutionMode(extract_mode);
  ExtractionDiagnostics diag;
  ExtractionDiagnostics* diag_ptr = diagnostics ? &diag : nullptr;

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
  R_xlen_t n = ids_zero_based.size();
  Rcpp::CharacterVector out_id(n);
  Rcpp::CharacterVector out_seq(n);
  Rcpp::CharacterVector out_qual;
  if (return_qual) {
    out_qual = Rcpp::CharacterVector(n);
  }

  std::unordered_map<std::uint64_t, std::vector<R_xlen_t>> positions_by_gid;
  positions_by_gid.reserve(static_cast<size_t>(requested.size()));
  for (R_xlen_t i = 0; i < n; ++i) {
    positions_by_gid[requested[static_cast<size_t>(i)]].push_back(i);
  }
  std::vector<std::uint64_t> unique_ids;
  unique_ids.reserve(positions_by_gid.size());
  for (const auto& kv : positions_by_gid) {
    unique_ids.push_back(kv.first);
  }
  std::sort(unique_ids.begin(), unique_ids.end());

  constexpr std::size_t kChunkUniqueIds = 50000;
  const bool need_parsed_qual = (bundle.format == "fastq" && include_qual);
  for (std::size_t start = 0; start < unique_ids.size(); start += kChunkUniqueIds) {
    const std::size_t end = std::min(start + kChunkUniqueIds, unique_ids.size());
    std::vector<std::uint64_t> chunk_ids(
      unique_ids.begin() + static_cast<std::ptrdiff_t>(start),
      unique_ids.begin() + static_cast<std::ptrdiff_t>(end)
    );
    const SelectedRecordMap selected_chunk = collectRequestedRecords(
      files,
      bundle,
      chunk_ids,
      need_parsed_qual,
      bridge_gap,
      region_bytes,
      mode,
      diag_ptr
    );
    for (std::uint64_t gid : chunk_ids) {
      auto rec_it = selected_chunk.find(gid);
      if (rec_it == selected_chunk.end()) {
        throw std::runtime_error("Could not resolve all requested ids with provided index/files.");
      }
      const ParsedRecord& rec = rec_it->second;
      SEXP id_ch = Rf_mkCharLen(rec.id.c_str(), static_cast<int>(rec.id.size()));
      SEXP seq_ch = Rf_mkCharLen(rec.seq.c_str(), static_cast<int>(rec.seq.size()));
      SEXP qual_ch = R_NilValue;
      if (return_qual) {
        qual_ch = Rf_mkCharLen(rec.qual.c_str(), static_cast<int>(rec.qual.size()));
      }

      const auto pos_it = positions_by_gid.find(gid);
      if (pos_it == positions_by_gid.end()) {
        continue;
      }
      for (R_xlen_t pos : pos_it->second) {
        SET_STRING_ELT(out_id, pos, id_ch);
        SET_STRING_ELT(out_seq, pos, seq_ch);
        if (return_qual) {
          SET_STRING_ELT(out_qual, pos, qual_ch);
        }
      }
    }
  }

  if (return_qual) {
    Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("seq_id") = out_id,
      Rcpp::Named("seq") = out_seq,
      Rcpp::Named("qual") = out_qual
    );
    if (diagnostics) {
      out.attr("fastqindexr_diagnostics") = diagnosticsToList(diag);
    }
    return out;
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("seq_id") = out_id,
    Rcpp::Named("seq") = out_seq
  );
  if (diagnostics) {
    out.attr("fastqindexr_diagnostics") = diagnosticsToList(diag);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_extract_sequences_to_file(
  Rcpp::CharacterVector files,
  std::string source_type,
  Rcpp::NumericVector ids_zero_based,
  SEXP index_ptr_sexp,
  std::string output_type,
  std::string outfile,
  bool append,
  bool compress,
  double max_bridge_gap,
  double max_region_bytes,
  std::string extract_mode,
  bool diagnostics
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
  const std::uint64_t bridge_gap = parseNonNegativeWhole(max_bridge_gap, "max_bridge_gap");
  const std::uint64_t region_bytes = parsePositiveWhole(max_region_bytes, "max_region_bytes");
  const ExtractExecutionMode mode = parseExtractExecutionMode(extract_mode);
  ExtractionDiagnostics diag;
  ExtractionDiagnostics* diag_ptr = diagnostics ? &diag : nullptr;
  if (ids_zero_based.size() < 1) {
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("written") = 0.0);
    if (diagnostics) {
      out.attr("fastqindexr_diagnostics") = diagnosticsToList(diag);
    }
    return out;
  }

  const std::string resolved_output_type = resolveOutputType(output_type, bundle.format);
  const std::vector<std::uint64_t> requested = parseRequestedIds(ids_zero_based);
  const bool need_parsed_qual = (resolved_output_type == "fastq");

  if (requestIsSortedUnique(requested) && !append) {
    const double written = streamSortedUniqueToFile(
      files,
      bundle,
      requested,
      resolved_output_type,
      outfile,
      compress,
      bridge_gap,
      region_bytes,
      mode,
      diag_ptr
    );
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("written") = written);
    if (diagnostics) {
      out.attr("fastqindexr_diagnostics") = diagnosticsToList(diag);
    }
    return out;
  }

  const SelectedRecordMap selected_global = collectRequestedRecords(
    files,
    bundle,
    requested,
    need_parsed_qual,
    bridge_gap,
    region_bytes,
    mode,
    diag_ptr
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
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("written") = static_cast<double>(written)
  );
  if (diagnostics) {
    out.attr("fastqindexr_diagnostics") = diagnosticsToList(diag);
  }
  return out;
}

// [[Rcpp::export]]
SEXP cpp_extract_sequences_dnastringset(
  Rcpp::CharacterVector files,
  std::string source_type,
  Rcpp::NumericVector ids_zero_based,
  SEXP index_ptr_sexp,
  double chunk_chars,
  double max_bridge_gap,
  double max_region_bytes,
  std::string extract_mode,
  bool diagnostics,
  std::string renumber_mode
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
  if (!std::isfinite(chunk_chars) || chunk_chars <= 0.0) {
    throw std::runtime_error("`chunk_chars` must be a finite positive value.");
  }
  const std::uint64_t bridge_gap = parseNonNegativeWhole(max_bridge_gap, "max_bridge_gap");
  const std::uint64_t region_bytes = parsePositiveWhole(max_region_bytes, "max_region_bytes");
  const ExtractExecutionMode mode = parseExtractExecutionMode(extract_mode);
  ExtractionDiagnostics diag;
  ExtractionDiagnostics* diag_ptr = diagnostics ? &diag : nullptr;
  if (renumber_mode != "none" &&
      renumber_mode != "zero_based" &&
      renumber_mode != "one_based") {
    throw std::runtime_error("Invalid renumber mode.");
  }
  if (ids_zero_based.size() == 0) {
    Rcpp::Environment biostrings = Rcpp::Environment::namespace_env("Biostrings");
    Rcpp::Function dna_string_set = biostrings["DNAStringSet"];
    SEXP empty = dna_string_set(Rcpp::CharacterVector(), Rcpp::Named("use.names") = true);
    if (diagnostics) {
      Rf_setAttrib(empty, Rf_install("fastqindexr_diagnostics"), diagnosticsToList(diag));
    }
    return empty;
  }

  const std::vector<std::uint64_t> requested = parseRequestedIds(ids_zero_based);
  SEXP out = buildDNAStringSetFromRequested(
    files,
    bundle,
    requested,
    static_cast<std::uint64_t>(chunk_chars),
    renumber_mode,
    bridge_gap,
    region_bytes,
    mode,
    diag_ptr
  );
  if (diagnostics) {
    Rf_setAttrib(out, Rf_install("fastqindexr_diagnostics"), diagnosticsToList(diag));
  }
  return out;
}
