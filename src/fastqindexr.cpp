#include <Rcpp.h>
#include <zlib.h>

#include "fastqindex_core/process/extract/Extractor.h"
#include "fastqindex_core/process/extract/IndexReader.h"
#include "fastqindex_core/process/index/Indexer.h"
#include "fastqindex_core/process/io/FileSource.h"

#include <algorithm>
#ifdef FASTQINDEXR_TIMING
#include <chrono>
#endif
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

#ifdef FASTQINDEXR_TIMING
constexpr bool kFastqindexrTimingEnabled = true;
#else
constexpr bool kFastqindexrTimingEnabled = false;
#endif

constexpr int kGzBufferSize = 8192;

struct ParsedRecord {
  std::string id;
  std::string seq;
  std::string qual;
};

enum class CompressionKind {
  Gzip = 0,
  Plain = 1
};

struct IndexedFileBase {
  std::uint64_t record_count{0};
  virtual ~IndexedFileBase() = default;
  virtual CompressionKind kind() const = 0;
  virtual std::unique_ptr<IndexedFileBase> clone() const = 0;
};

struct GzipIndexedFile : public IndexedFileBase {
  fastqindex_core::FileIndex index;

  CompressionKind kind() const override {
    return CompressionKind::Gzip;
  }
  std::unique_ptr<IndexedFileBase> clone() const override {
    return std::unique_ptr<IndexedFileBase>(new GzipIndexedFile(*this));
  }
};

struct PlainIndexedFile : public IndexedFileBase {
  std::vector<std::uint64_t> plain_header_byte_offsets;

  CompressionKind kind() const override {
    return CompressionKind::Plain;
  }
  std::unique_ptr<IndexedFileBase> clone() const override {
    return std::unique_ptr<IndexedFileBase>(new PlainIndexedFile(*this));
  }
};

struct IndexBundle {
  std::string format;
  int record_size{0};
  CompressionKind compression_kind{CompressionKind::Gzip};
  std::vector<std::string> files;
  std::vector<std::unique_ptr<IndexedFileBase>> indexed_files;
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

#ifdef FASTQINDEXR_TIMING
struct ExtractionDiagnostics {
  std::uint64_t regions_planned{0};
  std::uint64_t extract_attempts{0};
  std::uint64_t extract_failures{0};
  std::uint64_t extract_failures_partial{0};
  std::uint64_t extract_failures_malformed{0};
  std::uint64_t extract_failures_other{0};
  std::uint64_t fallback_invocations{0};
  std::uint64_t fallback_records{0};
  double time_region_plan_ms{0.0};
  double time_indexed_dense_ms{0.0};
  double time_indexed_selective_ms{0.0};
  double time_selective_seek_init_ms{0.0};
  double time_selective_inflate_ms{0.0};
  double time_selective_parse_ms{0.0};
  double time_selective_line_split_ms{0.0};
  double time_selective_materialize_ms{0.0};
  double time_selective_callback_ms{0.0};
  double time_fallback_ms{0.0};
  double time_sequential_init_ms{0.0};
  double time_sequential_inflate_ms{0.0};
  double time_sequential_parse_ms{0.0};
  double time_sequential_line_split_ms{0.0};
  double time_sequential_materialize_ms{0.0};
  double time_sequential_callback_ms{0.0};
  std::string last_failure_message{};
};

ExtractionDiagnostics* diagnosticsPtrIfEnabled(
  ExtractionDiagnostics* diag,
  bool diagnostics_requested
) {
  if constexpr (!kFastqindexrTimingEnabled) {
    return nullptr;
  }
  if (!diagnostics_requested) {
    return nullptr;
  }
  return diag;
}

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
    Rcpp::Named("time_region_plan_ms") = diag.time_region_plan_ms,
    Rcpp::Named("time_indexed_dense_ms") = diag.time_indexed_dense_ms,
    Rcpp::Named("time_indexed_selective_ms") = diag.time_indexed_selective_ms,
    Rcpp::Named("time_selective_seek_init_ms") = diag.time_selective_seek_init_ms,
    Rcpp::Named("time_selective_inflate_ms") = diag.time_selective_inflate_ms,
    Rcpp::Named("time_selective_parse_ms") = diag.time_selective_parse_ms,
    Rcpp::Named("time_selective_line_split_ms") = diag.time_selective_line_split_ms,
    Rcpp::Named("time_selective_materialize_ms") = diag.time_selective_materialize_ms,
    Rcpp::Named("time_selective_callback_ms") = diag.time_selective_callback_ms,
    Rcpp::Named("time_fallback_ms") = diag.time_fallback_ms,
    Rcpp::Named("time_sequential_init_ms") = diag.time_sequential_init_ms,
    Rcpp::Named("time_sequential_inflate_ms") = diag.time_sequential_inflate_ms,
    Rcpp::Named("time_sequential_parse_ms") = diag.time_sequential_parse_ms,
    Rcpp::Named("time_sequential_line_split_ms") = diag.time_sequential_line_split_ms,
    Rcpp::Named("time_sequential_materialize_ms") = diag.time_sequential_materialize_ms,
    Rcpp::Named("time_sequential_callback_ms") = diag.time_sequential_callback_ms,
    Rcpp::Named("last_failure_message") = diag.last_failure_message
  );
}

#endif

std::uint64_t asUInt64(double value, const std::string& field) {
  if (!std::isfinite(value) || value < 0) {
    throw std::runtime_error("Malformed serialized index field: " + field);
  }
  return static_cast<std::uint64_t>(value);
}

const GzipIndexedFile& asGzipIndexedFile(const IndexedFileBase& indexed_file) {
  const auto* ptr = dynamic_cast<const GzipIndexedFile*>(&indexed_file);
  if (ptr == nullptr) {
    throw std::runtime_error("Expected gzip indexed file.");
  }
  return *ptr;
}

const PlainIndexedFile& asPlainIndexedFile(const IndexedFileBase& indexed_file) {
  const auto* ptr = dynamic_cast<const PlainIndexedFile*>(&indexed_file);
  if (ptr == nullptr) {
    throw std::runtime_error("Expected plain indexed file.");
  }
  return *ptr;
}

bool isPlainIndexedFile(const IndexedFileBase& indexed_file) {
  return indexed_file.kind() == CompressionKind::Plain;
}

Rcpp::List serializeIndexedFile(const IndexedFileBase& indexed_file) {
  const bool plain = isPlainIndexedFile(indexed_file);
  const auto* gzip = plain ? nullptr : &asGzipIndexedFile(indexed_file);
  const auto* plain_file = plain ? &asPlainIndexedFile(indexed_file) : nullptr;
  const auto& entries = plain ? std::vector<fastqindex_core::IndexEntry>() : gzip->index.entries;
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

  Rcpp::NumericVector plain_off(
    plain_file == nullptr ? 0 : static_cast<R_xlen_t>(plain_file->plain_header_byte_offsets.size())
  );
  for (R_xlen_t i = 0; i < plain_off.size(); ++i) {
    plain_off[i] = static_cast<double>(plain_file->plain_header_byte_offsets[static_cast<size_t>(i)]);
  }

  return Rcpp::List::create(
    Rcpp::Named("block_index") = block_index,
    Rcpp::Named("block_offset_in_raw_file") = block_offset,
    Rcpp::Named("starting_line_in_entry") = starting_line,
    Rcpp::Named("offset_to_next_line_start") = offset_next_line,
    Rcpp::Named("compressed_dictionary_size") = compressed_dictionary_size,
    Rcpp::Named("bits") = bits,
    Rcpp::Named("dictionary_blob") = dictionary_blob,
    Rcpp::Named("total_lines") = plain ? 0.0 :
      static_cast<double>(gzip->index.total_lines),
    Rcpp::Named("record_count") = static_cast<double>(indexed_file.record_count),
    Rcpp::Named("plain") = plain,
    Rcpp::Named("plain_header_byte_offsets") = plain_off
  );
}

Rcpp::List serializeIndexBundle(const IndexBundle& bundle) {
  Rcpp::List indexed_files(bundle.indexed_files.size());
  for (R_xlen_t i = 0; i < indexed_files.size(); ++i) {
    indexed_files[i] = serializeIndexedFile(*bundle.indexed_files[static_cast<size_t>(i)]);
  }

  Rcpp::NumericVector record_offsets(bundle.record_offsets.size());
  for (R_xlen_t i = 0; i < record_offsets.size(); ++i) {
    record_offsets[i] = static_cast<double>(bundle.record_offsets[static_cast<size_t>(i)]);
  }

  const int schema_out = 4;
  return Rcpp::List::create(
    Rcpp::Named("schema_version") = schema_out,
    Rcpp::Named("format") = bundle.format,
    Rcpp::Named("record_size") = bundle.record_size,
    Rcpp::Named("bundle_compression") = (bundle.compression_kind == CompressionKind::Gzip) ? "gzip" : "plain",
    Rcpp::Named("files") = bundle.files,
    Rcpp::Named("record_offsets") = record_offsets,
    Rcpp::Named("n_records") = static_cast<double>(bundle.n_records),
    Rcpp::Named("indexed_files") = indexed_files
  );
}

std::unique_ptr<IndexedFileBase> deserializeIndexedFile(
  const Rcpp::List& payload,
  int schema_version,
  CompressionKind bundle_kind,
  int record_size
) {
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

  const bool payload_plain = payload.containsElementNamed("plain") && Rcpp::as<bool>(payload["plain"]);
  if (
    (bundle_kind == CompressionKind::Gzip && payload_plain) ||
      (bundle_kind == CompressionKind::Plain && !payload_plain)
  ) {
    throw std::runtime_error("Malformed payload: indexed file compression differs from bundle compression.");
  }
  if (bundle_kind == CompressionKind::Gzip) {
    std::unique_ptr<GzipIndexedFile> indexed_file(new GzipIndexedFile());
    indexed_file->index.entries.reserve(static_cast<size_t>(n));
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
      indexed_file->index.entries.push_back(std::move(entry));
    }
    indexed_file->index.total_lines = asUInt64(Rcpp::as<double>(payload["total_lines"]), "total_lines");
    indexed_file->record_count = asUInt64(Rcpp::as<double>(payload["record_count"]), "record_count");
    return indexed_file;
  }

  std::unique_ptr<PlainIndexedFile> indexed_file(new PlainIndexedFile());
  if (!payload.containsElementNamed("plain_header_byte_offsets")) {
    throw std::runtime_error("Malformed payload: plain indexed file missing plain_header_byte_offsets.");
  }
  Rcpp::NumericVector po = payload["plain_header_byte_offsets"];
  indexed_file->plain_header_byte_offsets.reserve(static_cast<size_t>(po.size()));
  for (R_xlen_t i = 0; i < po.size(); ++i) {
    indexed_file->plain_header_byte_offsets.push_back(asUInt64(po[i], "plain_header_byte_offsets"));
  }
  indexed_file->record_count = asUInt64(Rcpp::as<double>(payload["record_count"]), "record_count");
  if (static_cast<int>(indexed_file->plain_header_byte_offsets.size()) != static_cast<int>(indexed_file->record_count)) {
    throw std::runtime_error("Malformed payload: plain_header_byte_offsets length must equal record_count.");
  }
  return indexed_file;
}

IndexBundle deserializeIndexBundle(const Rcpp::List& payload) {
  const int schema_version = Rcpp::as<int>(payload["schema_version"]);
  if (schema_version != 4) {
    throw std::runtime_error(
      "Unsupported serialized index schema version. This refactor requires rebuilding indexes."
    );
  }

  IndexBundle bundle;
  bundle.format = Rcpp::as<std::string>(payload["format"]);
  bundle.record_size = Rcpp::as<int>(payload["record_size"]);
  std::string bundle_compression = Rcpp::as<std::string>(payload["bundle_compression"]);
  if (bundle_compression == "gzip") {
    bundle.compression_kind = CompressionKind::Gzip;
  } else if (bundle_compression == "plain") {
    bundle.compression_kind = CompressionKind::Plain;
  } else {
    throw std::runtime_error("Malformed payload: unsupported bundle_compression.");
  }
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
    bundle.indexed_files.push_back(
      deserializeIndexedFile(indexed_files[i], schema_version, bundle.compression_kind, bundle.record_size)
    );
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

bool fileHasGzipMagic(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in.good()) {
    return false;
  }
  unsigned char b0 = 0;
  unsigned char b1 = 0;
  in.read(reinterpret_cast<char*>(&b0), 1);
  in.read(reinterpret_cast<char*>(&b1), 1);
  return in.good() && b0 == 0x1f && b1 == 0x8b;
}

std::string trimPrefix(const std::string& line, char prefix) {
  if (!line.empty() && line.front() == prefix) {
    return line.substr(1);
  }
  return line;
}

std::string openAndDetectType(const std::string& path) {
  if (fileHasGzipMagic(path)) {
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
  std::ifstream in(path, std::ios::binary);
  if (!in.good()) {
    throw std::runtime_error("Could not open input: " + path);
  }
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    if (line.empty()) {
      continue;
    }
    if (line.front() == '>') {
      return "fasta";
    }
    if (line.front() == '@') {
      return "fastq";
    }
    throw std::runtime_error("Cannot infer type from first non-empty line in: " + path);
  }
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

std::unique_ptr<GzipIndexedFile> readSingleFqi(const std::string& fqi_file) {
  auto source = fastqindex_core::FileSource::from(std::filesystem::path(fqi_file));
  fastqindex_core::IndexReader reader(source);
  std::vector<std::shared_ptr<fastqindex_core::IndexEntry>> upstream_entries =
    reader.readIndexFile();
  if (upstream_entries.empty()) {
    const auto errors = reader.getErrorMessages();
    std::string suffix = errors.empty() ? "" : (" (" + errors.front() + ")");
    throw std::runtime_error("Failed to parse .fqi index: " + fqi_file + suffix);
  }

  std::unique_ptr<GzipIndexedFile> indexed_file(new GzipIndexedFile());
  indexed_file->index.entries.reserve(upstream_entries.size());
  for (const auto& upstream_entry : upstream_entries) {
    indexed_file->index.entries.push_back(*upstream_entry);
  }
  indexed_file->index.total_lines = static_cast<std::uint64_t>(
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
  std::uint64_t max_region_records,
  int record_size
  #ifdef FASTQINDEXR_TIMING
  ,
  ExtractionDiagnostics* diag
  #endif
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
  for (size_t i = 1; i < ids.size(); ++i) {
    const bool within_gap = (ids[i] <= last + max_bridge_gap + 1);
    const std::uint64_t span_records = ids[i] - start + 1;
    const bool within_records = (span_records <= max_region_records);
    if (within_gap && within_records) {
      last = ids[i];
    } else {
      regions.emplace_back(start, last);
      start = ids[i];
      last = ids[i];
    }
  }
  regions.emplace_back(start, last);
  #ifdef FASTQINDEXR_TIMING
  if (diag != nullptr) {
    diag->regions_planned += static_cast<std::uint64_t>(regions.size());
  }
  #endif
  return regions;
}

std::vector<std::pair<std::uint64_t, std::uint64_t>> denseRegions(
  std::vector<std::uint64_t> ids,
  std::uint64_t max_bridge_gap,
  std::uint64_t max_region_records,
  int record_size
  #ifdef FASTQINDEXR_TIMING
  ,
  ExtractionDiagnostics* diag
  #endif
) {
  if (ids.empty()) {
    return {};
  }
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
  return buildDenseRegionsFromSortedUnique(
    ids,
    max_bridge_gap,
    max_region_records,
      record_size
      #ifdef FASTQINDEXR_TIMING
      ,
      diag
      #endif
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

std::uint64_t scanSequentialPlainRequestedRecords(
  const std::string& file,
  const std::string& type,
  bool include_qual,
  const std::unordered_set<std::uint64_t>& requested_local_ids,
  const std::function<void(std::uint64_t, const ParsedRecord&)>& on_record
  #ifdef FASTQINDEXR_TIMING
  ,
  ExtractionDiagnostics* diag
  #endif
) {
  #ifdef FASTQINDEXR_TIMING
  if (diag != nullptr) {
    diag->fallback_invocations++;
  }
  const auto init_t0 = std::chrono::steady_clock::now();
  #endif
  std::ifstream stream(file, std::ios::binary);
  if (!stream.good()) {
    throw std::runtime_error("Fallback parser could not open plain file: " + file);
  }
  #ifdef FASTQINDEXR_TIMING
  if (diag != nullptr) {
    const auto init_t1 = std::chrono::steady_clock::now();
    diag->time_sequential_init_ms +=
      std::chrono::duration<double, std::milli>(init_t1 - init_t0).count();
  }
  #endif

  std::vector<std::uint64_t> requested_sorted(requested_local_ids.begin(), requested_local_ids.end());
  std::sort(requested_sorted.begin(), requested_sorted.end());
  std::size_t requested_idx = 0;
  auto is_requested = [&](std::uint64_t local_id) -> bool {
    while (requested_idx < requested_sorted.size() && requested_sorted[requested_idx] < local_id) {
      requested_idx++;
    }
    return requested_idx < requested_sorted.size() && requested_sorted[requested_idx] == local_id;
  };

  const bool is_fastq = (type == "fastq");
  const int lines_per_record = is_fastq ? 4 : 2;
  std::uint64_t local_id = 0;
  std::uint64_t selected_count = 0;
  bool capture_current_record = is_requested(local_id);
  int line_in_record = 0;
  std::string current_line;
  current_line.reserve(256);
  std::string header;
  std::string seq;
  std::string qual;
  std::vector<char> buffer(kGzBufferSize * 8, '\0');
  const std::string empty_line;

  auto consume_line = [&](const std::string& line) {
    #ifdef FASTQINDEXR_TIMING
    const auto parse_t0 = std::chrono::steady_clock::now();
    #endif
    if (capture_current_record) {
      if (is_fastq) {
        if (line_in_record == 0) {
          header = line;
        } else if (line_in_record == 1) {
          seq = line;
        } else if (line_in_record == 3) {
          qual = line;
        }
      } else {
        if (line_in_record == 0) {
          header = line;
        } else if (line_in_record == 1) {
          seq = line;
        }
      }
    }
    line_in_record++;
    if (line_in_record == lines_per_record) {
      #ifdef FASTQINDEXR_TIMING
      if (diag != nullptr) {
        const auto parse_t1 = std::chrono::steady_clock::now();
        diag->time_sequential_parse_ms +=
          std::chrono::duration<double, std::milli>(parse_t1 - parse_t0).count();
      }
      #endif
      if (capture_current_record) {
        #ifdef FASTQINDEXR_TIMING
        const auto materialize_t0 = std::chrono::steady_clock::now();
        #endif
        const ParsedRecord rec{
          trimPrefix(header, is_fastq ? '@' : '>'),
          seq,
          (is_fastq && include_qual) ? qual : ""
        };
        #ifdef FASTQINDEXR_TIMING
        if (diag != nullptr) {
          const auto materialize_t1 = std::chrono::steady_clock::now();
          diag->time_sequential_materialize_ms +=
            std::chrono::duration<double, std::milli>(materialize_t1 - materialize_t0).count();
        }
        const auto callback_t0 = std::chrono::steady_clock::now();
        #endif
        on_record(local_id, rec);
        #ifdef FASTQINDEXR_TIMING
        if (diag != nullptr) {
          const auto callback_t1 = std::chrono::steady_clock::now();
          diag->time_sequential_callback_ms +=
            std::chrono::duration<double, std::milli>(callback_t1 - callback_t0).count();
          diag->fallback_records++;
        }
        #endif
        selected_count++;
      }
      local_id++;
      capture_current_record = is_requested(local_id);
      line_in_record = 0;
      header.clear();
      seq.clear();
      qual.clear();
      return;
    }
    #ifdef FASTQINDEXR_TIMING
    if (diag != nullptr) {
      const auto parse_t1 = std::chrono::steady_clock::now();
      diag->time_sequential_parse_ms +=
        std::chrono::duration<double, std::milli>(parse_t1 - parse_t0).count();
    }
    #endif
  };

  while (true) {
    stream.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    const std::streamsize n_read_s = stream.gcount();
    if (n_read_s <= 0) {
      break;
    }
    const auto n_read = static_cast<std::size_t>(n_read_s);
    #ifdef FASTQINDEXR_TIMING
    const auto split_t0 = std::chrono::steady_clock::now();
    #endif
    std::size_t line_start = 0;
    for (std::size_t i = 0; i < n_read; ++i) {
      const char ch = buffer[i];
      if (ch != '\n') {
        continue;
      }
      std::size_t seg_end = i;
      if (seg_end > line_start && buffer[seg_end - 1] == '\r') {
        seg_end--;
      }
      if (capture_current_record) {
        if (seg_end > line_start) {
          const char* seg_ptr = buffer.data() + line_start;
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
        consume_line(empty_line);
      }
      line_start = i + 1;
    }
    if (line_start < n_read && capture_current_record) {
      const char* seg_ptr = buffer.data() + line_start;
      const std::size_t seg_len = n_read - line_start;
      if (current_line.empty()) {
        current_line.assign(seg_ptr, seg_len);
      } else {
        current_line.append(seg_ptr, seg_len);
      }
    }
    #ifdef FASTQINDEXR_TIMING
    if (diag != nullptr) {
      const auto split_t1 = std::chrono::steady_clock::now();
      diag->time_sequential_line_split_ms +=
        std::chrono::duration<double, std::milli>(split_t1 - split_t0).count();
    }
    #endif
  }
  if (!current_line.empty()) {
    consume_line(current_line);
  }

  return selected_count;
}

std::uint64_t scanSequentialRequestedRecords(
  const std::string& file,
  const std::string& type,
  bool include_qual,
  const std::unordered_set<std::uint64_t>& requested_local_ids,
  const std::function<void(std::uint64_t, const ParsedRecord&)>& on_record
  #ifdef FASTQINDEXR_TIMING
  ,
  ExtractionDiagnostics* diag
  #endif
) {
  if (!fileHasGzipMagic(file)) {
    return scanSequentialPlainRequestedRecords(
      file,
      type,
      include_qual,
      requested_local_ids,
      on_record
      #ifdef FASTQINDEXR_TIMING
      ,
      diag
      #endif
    );
  }
  #ifdef FASTQINDEXR_TIMING
  if (diag != nullptr) {
    diag->fallback_invocations++;
  }
  const auto init_t0 = std::chrono::steady_clock::now();
  #endif
  gzFile stream = gzopen(file.c_str(), "rb");
  if (stream == nullptr) {
    throw std::runtime_error("Fallback parser could not open file: " + file);
  }
  #ifdef FASTQINDEXR_TIMING
  if (diag != nullptr) {
    const auto init_t1 = std::chrono::steady_clock::now();
    diag->time_sequential_init_ms +=
      std::chrono::duration<double, std::milli>(init_t1 - init_t0).count();
  }
  #endif

  std::vector<std::uint64_t> requested_sorted(requested_local_ids.begin(), requested_local_ids.end());
  std::sort(requested_sorted.begin(), requested_sorted.end());
  std::size_t requested_idx = 0;
  auto is_requested = [&](std::uint64_t local_id) -> bool {
    while (requested_idx < requested_sorted.size() && requested_sorted[requested_idx] < local_id) {
      requested_idx++;
    }
    return requested_idx < requested_sorted.size() && requested_sorted[requested_idx] == local_id;
  };

  const bool is_fastq = (type == "fastq");
  const int lines_per_record = is_fastq ? 4 : 2;
  std::uint64_t local_id = 0;
  std::uint64_t selected_count = 0;
  bool capture_current_record = is_requested(local_id);
  int line_in_record = 0;
  std::string current_line;
  current_line.reserve(256);
  std::string header;
  std::string seq;
  std::string qual;
  std::vector<char> buffer(kGzBufferSize * 8, '\0');
  const std::string empty_line;

  auto consume_line = [&](const std::string& line) {
    #ifdef FASTQINDEXR_TIMING
    const auto parse_t0 = std::chrono::steady_clock::now();
    #endif
    if (capture_current_record) {
      if (is_fastq) {
        if (line_in_record == 0) {
          header = line;
        } else if (line_in_record == 1) {
          seq = line;
        } else if (line_in_record == 3) {
          qual = line;
        }
      } else {
        if (line_in_record == 0) {
          header = line;
        } else if (line_in_record == 1) {
          seq = line;
        }
      }
    }
    line_in_record++;
    if (line_in_record == lines_per_record) {
      #ifdef FASTQINDEXR_TIMING
      if (diag != nullptr) {
        const auto parse_t1 = std::chrono::steady_clock::now();
        diag->time_sequential_parse_ms +=
          std::chrono::duration<double, std::milli>(parse_t1 - parse_t0).count();
      }
      #endif
      if (capture_current_record) {
        #ifdef FASTQINDEXR_TIMING
        const auto materialize_t0 = std::chrono::steady_clock::now();
        #endif
        const ParsedRecord rec{
          trimPrefix(header, is_fastq ? '@' : '>'),
          seq,
          (is_fastq && include_qual) ? qual : ""
        };
        #ifdef FASTQINDEXR_TIMING
        if (diag != nullptr) {
          const auto materialize_t1 = std::chrono::steady_clock::now();
          diag->time_sequential_materialize_ms +=
            std::chrono::duration<double, std::milli>(materialize_t1 - materialize_t0).count();
        }
        const auto callback_t0 = std::chrono::steady_clock::now();
        #endif
        on_record(local_id, rec);
        #ifdef FASTQINDEXR_TIMING
        if (diag != nullptr) {
          const auto callback_t1 = std::chrono::steady_clock::now();
          diag->time_sequential_callback_ms +=
            std::chrono::duration<double, std::milli>(callback_t1 - callback_t0).count();
          diag->fallback_records++;
        }
        #endif
        selected_count++;
      }
      local_id++;
      capture_current_record = is_requested(local_id);
      line_in_record = 0;
      header.clear();
      seq.clear();
      qual.clear();
      return;
    }
    #ifdef FASTQINDEXR_TIMING
    if (diag != nullptr) {
      const auto parse_t1 = std::chrono::steady_clock::now();
      diag->time_sequential_parse_ms +=
        std::chrono::duration<double, std::milli>(parse_t1 - parse_t0).count();
    }
    #endif
  };

  while (true) {
    #ifdef FASTQINDEXR_TIMING
    const auto inflate_t0 = std::chrono::steady_clock::now();
    #endif
    const int n_read = gzread(stream, buffer.data(), static_cast<unsigned int>(buffer.size()));
    #ifdef FASTQINDEXR_TIMING
    if (diag != nullptr) {
      const auto inflate_t1 = std::chrono::steady_clock::now();
      diag->time_sequential_inflate_ms +=
        std::chrono::duration<double, std::milli>(inflate_t1 - inflate_t0).count();
    }
    #endif
    if (n_read <= 0) {
      break;
    }
    #ifdef FASTQINDEXR_TIMING
    const auto split_t0 = std::chrono::steady_clock::now();
    #endif
    std::size_t line_start = 0;
    for (std::size_t i = 0; i < static_cast<std::size_t>(n_read); ++i) {
      const char ch = buffer[i];
      if (ch != '\n') {
        continue;
      }
      std::size_t seg_end = i;
      if (seg_end > line_start && buffer[seg_end - 1] == '\r') {
        seg_end--;
      }
      if (capture_current_record) {
        if (seg_end > line_start) {
          const char* seg_ptr = buffer.data() + line_start;
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
        consume_line(empty_line);
      }
      line_start = i + 1;
    }
    if (line_start < static_cast<std::size_t>(n_read) && capture_current_record) {
      const char* seg_ptr = buffer.data() + line_start;
      const std::size_t seg_len = static_cast<std::size_t>(n_read) - line_start;
      if (current_line.empty()) {
        current_line.assign(seg_ptr, seg_len);
      } else {
        current_line.append(seg_ptr, seg_len);
      }
    }
    #ifdef FASTQINDEXR_TIMING
    if (diag != nullptr) {
      const auto split_t1 = std::chrono::steady_clock::now();
      diag->time_sequential_line_split_ms +=
        std::chrono::duration<double, std::milli>(split_t1 - split_t0).count();
    }
    #endif
  }
  if (!current_line.empty()) {
    consume_line(current_line);
  }

  gzclose(stream);
  return selected_count;
}

void sequentialFallbackCollect(
  const std::string& file,
  const std::string& type,
  const std::unordered_set<std::uint64_t>& requested_local_ids,
  std::uint64_t global_offset,
  std::unordered_map<std::uint64_t, ParsedRecord>* selected_global,
  bool include_qual
  #ifdef FASTQINDEXR_TIMING
  , ExtractionDiagnostics* diag
  #endif
) {
  #ifdef FASTQINDEXR_TIMING
  const auto t0 = std::chrono::steady_clock::now();
  #endif
  scanSequentialRequestedRecords(
    file,
    type,
    include_qual,
    requested_local_ids,
    [&](std::uint64_t local_id, const ParsedRecord& rec) {
      (*selected_global)[global_offset + local_id] = rec;
    }
    #ifdef FASTQINDEXR_TIMING
    , diag
    #endif
  );
  #ifdef FASTQINDEXR_TIMING
  if (diag != nullptr) {
    const auto t1 = std::chrono::steady_clock::now();
    diag->time_fallback_ms +=
      std::chrono::duration<double, std::milli>(t1 - t0).count();
  }
  #endif
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
  if (mode == "sequential_only" || mode == "sequential") {
    return ExtractExecutionMode::SequentialOnly;
  }
  throw std::runtime_error("`extract_mode` must be 'indexed', 'sequential', or 'sequential_only'.");
}

ParsedRecord readOneRecordFromPlainStream(
  std::istream* in,
  const std::string& type,
  bool include_qual
) {
  std::string line;
  std::string header;
  std::string seq;
  std::string qual;
  const bool is_fastq = (type == "fastq");
  const int lines_per_record = is_fastq ? 4 : 2;
  for (int li = 0; li < lines_per_record; ++li) {
    if (!std::getline(*in, line)) {
      throw std::runtime_error("Unexpected EOF while reading plain FASTA/FASTQ record.");
    }
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    if (li == 0) {
      header = line;
    } else if (li == 1) {
      seq = line;
    } else if (is_fastq && li == 3) {
      qual = line;
    }
  }
  const char prefix = is_fastq ? '@' : '>';
  return ParsedRecord{
    trimPrefix(header, prefix),
    seq,
    (is_fastq && include_qual) ? qual : ""
  };
}

ParsedRecord readPlainRecordAtOffset(
  const std::string& path,
  std::uint64_t byte_offset,
  const std::string& type,
  bool include_qual
) {
  std::ifstream in(path, std::ios::binary);
  if (!in.good()) {
    throw std::runtime_error("Could not open plain input: " + path);
  }
  in.seekg(static_cast<std::streamoff>(byte_offset), std::ios::beg);
  if (!in.good()) {
    throw std::runtime_error("Could not seek plain input: " + path);
  }
  return readOneRecordFromPlainStream(&in, type, include_qual);
}

std::unique_ptr<PlainIndexedFile> buildPlainFileIndex(const std::string& path, int record_size) {
  std::unique_ptr<PlainIndexedFile> out(new PlainIndexedFile());
  std::ifstream in(path, std::ios::binary);
  if (!in.good()) {
    throw std::runtime_error("Could not open plain input for indexing: " + path);
  }
  std::uint64_t line_index = 0;
  std::string line;
  while (true) {
    const std::streampos line_start = in.tellg();
    if (!std::getline(in, line)) {
      break;
    }
    if (line_index % static_cast<std::uint64_t>(record_size) == 0) {
      if (line_start < 0) {
        throw std::runtime_error("Could not track byte position in plain file: " + path);
      }
      out->plain_header_byte_offsets.push_back(static_cast<std::uint64_t>(line_start));
    }
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    line_index++;
  }
  if (line_index % static_cast<std::uint64_t>(record_size) != 0) {
    throw std::runtime_error("Plain file line count is not a multiple of record size for: " + path);
  }
  out->record_count = line_index / static_cast<std::uint64_t>(record_size);
  if (out->plain_header_byte_offsets.size() != out->record_count) {
    throw std::runtime_error("Plain index header offset count mismatch for: " + path);
  }
  return out;
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
  std::uint64_t max_region_records,
  ExtractExecutionMode mode
  #ifdef FASTQINDEXR_TIMING
  , ExtractionDiagnostics* diag
  #endif
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
        need_qual
        #ifdef FASTQINDEXR_TIMING
        , diag
        #endif
      );
      continue;
    }
    if (isPlainIndexedFile(*bundle.indexed_files[file_idx]) && mode == ExtractExecutionMode::Indexed) {
      const auto& offs = asPlainIndexedFile(*bundle.indexed_files[file_idx]).plain_header_byte_offsets;
      for (auto lid : local_ids) {
        if (lid >= offs.size()) {
          throw std::runtime_error("Requested id out of range for plain indexed file.");
        }
        ParsedRecord rec = readPlainRecordAtOffset(
          Rcpp::as<std::string>(files[file_idx]),
          offs[static_cast<size_t>(lid)],
          bundle.format,
          need_qual
        );
        selected_global[bundle.record_offsets[file_idx] + lid] = std::move(rec);
      }
      continue;
    }
    #ifdef FASTQINDEXR_TIMING
    const auto plan_t0 = std::chrono::steady_clock::now();
    #endif
    auto regions = denseRegions(
      local_ids,
      max_bridge_gap,
      max_region_records,
      bundle.record_size
      #ifdef FASTQINDEXR_TIMING
      , diag
      #endif
    );
    #ifdef FASTQINDEXR_TIMING
    if (diag != nullptr) {
      const auto plan_t1 = std::chrono::steady_clock::now();
      diag->time_region_plan_ms +=
        std::chrono::duration<double, std::milli>(plan_t1 - plan_t0).count();
    }
    #endif
    for (const auto& rg : regions) {
      std::uint64_t region_start = rg.first;
      std::uint64_t region_end = rg.second;
      std::uint64_t line_start = region_start * static_cast<std::uint64_t>(bundle.record_size);
      std::uint64_t line_count =
        (region_end - region_start + 1) * static_cast<std::uint64_t>(bundle.record_size);

      try {
        #ifdef FASTQINDEXR_TIMING
        if (diag != nullptr) {
          diag->extract_attempts++;
        }
        #endif
        const auto region_begin = std::lower_bound(local_ids.begin(), local_ids.end(), region_start);
        const auto region_end_it = std::upper_bound(local_ids.begin(), local_ids.end(), region_end);
        const std::uint64_t requested_in_region = static_cast<std::uint64_t>(region_end_it - region_begin);
        const std::uint64_t region_span = region_end - region_start + 1;
        const bool dense_region = (requested_in_region == region_span);

        if (dense_region) {
          #ifdef FASTQINDEXR_TIMING
          const auto dense_t0 = std::chrono::steady_clock::now();
          #endif
          std::uint64_t extracted_records = 0;
          extractor.extractRecords(
            Rcpp::as<std::string>(files[file_idx]),
            asGzipIndexedFile(*bundle.indexed_files[file_idx]).index.entries,
            line_start,
            line_count,
            bundle.record_size,
            region_start,
            bundle.format,
            need_qual,
            [&](std::uint64_t local_id, const fastqindex_core::ExtractedRecord& rec_in) {
              ParsedRecord rec{rec_in.seq_id, rec_in.seq, rec_in.qual};
              const std::uint64_t global_id = bundle.record_offsets[file_idx] + local_id;
              selected_global[global_id] = std::move(rec);
              extracted_records++;
            }
          );
          if (extracted_records < region_span) {
            throw std::runtime_error("Partial extraction for requested region.");
          }
          if (extracted_records > region_span) {
            throw std::runtime_error("Malformed extracted line count for region.");
          }
          #ifdef FASTQINDEXR_TIMING
          if (diag != nullptr) {
            const auto dense_t1 = std::chrono::steady_clock::now();
            diag->time_indexed_dense_ms +=
              std::chrono::duration<double, std::milli>(dense_t1 - dense_t0).count();
          }
          #endif
        } else {
          #ifdef FASTQINDEXR_TIMING
          const auto sel_t0 = std::chrono::steady_clock::now();
          #endif
          const std::unordered_set<std::uint64_t> region_requested_set(region_begin, region_end_it);
          #ifdef FASTQINDEXR_TIMING
          fastqindex_core::SelectiveExtractTimings selective_timings;
          fastqindex_core::SelectiveExtractTimings* selective_timings_ptr =
            (diag != nullptr) ? &selective_timings : nullptr;
          #endif
          std::uint64_t selected_records = 0;
          extractor.extractSelectedRecords(
            Rcpp::as<std::string>(files[file_idx]),
            asGzipIndexedFile(*bundle.indexed_files[file_idx]).index.entries,
            line_start,
            line_count,
            bundle.record_size,
            region_requested_set,
            region_start,
            bundle.format,
            need_qual,
            [&](std::uint64_t local_id, const fastqindex_core::ExtractedRecord& rec_in) {
              ParsedRecord rec{rec_in.seq_id, rec_in.seq, rec_in.qual};
              const std::uint64_t global_id = bundle.record_offsets[file_idx] + local_id;
              selected_global[global_id] = std::move(rec);
              selected_records++;
            }
            #ifdef FASTQINDEXR_TIMING
            ,
            selective_timings_ptr
            #endif
          );
          if (selected_records < requested_in_region) {
            throw std::runtime_error("Partial extraction for requested region.");
          }
          #ifdef FASTQINDEXR_TIMING
          if (diag != nullptr) {
            const auto sel_t1 = std::chrono::steady_clock::now();
            diag->time_indexed_selective_ms +=
              std::chrono::duration<double, std::milli>(sel_t1 - sel_t0).count();
            diag->time_selective_seek_init_ms += selective_timings_ptr->time_seek_init_ms;
            diag->time_selective_inflate_ms += selective_timings_ptr->time_inflate_ms;
            diag->time_selective_parse_ms += selective_timings_ptr->time_parse_ms;
            diag->time_selective_line_split_ms += selective_timings_ptr->time_line_split_ms;
            diag->time_selective_materialize_ms += selective_timings_ptr->time_materialize_ms;
            diag->time_selective_callback_ms += selective_timings_ptr->time_callback_ms;
          }
          #endif
        }
      } catch (const std::exception& e) {
        #ifdef FASTQINDEXR_TIMING
        if (diag != nullptr) {
          recordExtractFailure(diag, e.what());
        }
        #endif
        // fastqindexr change: Fallback preserves correctness when random-seek
        // extraction fails for particular gzip layouts.
        sequentialFallbackCollect(
          Rcpp::as<std::string>(files[file_idx]),
          bundle.format,
          requested_set,
          bundle.record_offsets[file_idx],
          &selected_global,
          need_qual
          #ifdef FASTQINDEXR_TIMING
          , diag
          #endif
        );
        break;
      } catch (...) {
        #ifdef FASTQINDEXR_TIMING
        recordExtractFailure(diag, "Unknown indexed extraction failure.");
        #endif
        sequentialFallbackCollect(
          Rcpp::as<std::string>(files[file_idx]),
          bundle.format,
          requested_set,
          bundle.record_offsets[file_idx],
          &selected_global,
          need_qual
          #ifdef FASTQINDEXR_TIMING
          , diag
          #endif
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
  OutputWriter* writer
  #ifdef FASTQINDEXR_TIMING
  , ExtractionDiagnostics* diag
  #endif
) {
  #ifdef FASTQINDEXR_TIMING
  const auto t0 = std::chrono::steady_clock::now();
  #endif
  const std::uint64_t count = scanSequentialRequestedRecords(
    file,
    type,
    need_qual,
    requested_local_ids,
    [&](std::uint64_t, const ParsedRecord& rec) {
      writer->write(renderRecord(rec, output_type));
    }
    #ifdef FASTQINDEXR_TIMING
    , diag
    #endif
  );
  #ifdef FASTQINDEXR_TIMING
  if (diag != nullptr) {
    const auto t1 = std::chrono::steady_clock::now();
    diag->time_fallback_ms +=
      std::chrono::duration<double, std::milli>(t1 - t0).count();
  }
  #endif
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
  std::uint64_t max_region_records,
  ExtractExecutionMode mode
  #ifdef FASTQINDEXR_TIMING
  , ExtractionDiagnostics* diag
  #endif
) {
  const bool include_qual = (resolved_output_type == "fastq");
  const SelectedRecordMap selected_global = collectRequestedRecords(
    files,
    bundle,
    requested,
    include_qual,
    max_bridge_gap,
    max_region_records,
    mode
    #ifdef FASTQINDEXR_TIMING
    , diag
    #endif
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
  std::uint64_t max_region_records,
  ExtractExecutionMode mode
  #ifdef FASTQINDEXR_TIMING
  , ExtractionDiagnostics* diag
  #endif
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
        writer.get()
        #ifdef FASTQINDEXR_TIMING
        , diag
        #endif
      );
      continue;
    }
    if (isPlainIndexedFile(*bundle.indexed_files[file_idx]) && mode == ExtractExecutionMode::Indexed) {
      const auto& offs = asPlainIndexedFile(*bundle.indexed_files[file_idx]).plain_header_byte_offsets;
      for (auto lid : local_ids) {
        if (lid >= offs.size()) {
          throw std::runtime_error("Requested id out of range for plain indexed file.");
        }
        ParsedRecord rec = readPlainRecordAtOffset(
          Rcpp::as<std::string>(files[file_idx]),
          offs[static_cast<size_t>(lid)],
          bundle.format,
          need_qual
        );
        writer->write(renderRecord(rec, resolved_output_type));
        n_written++;
      }
      continue;
    }
    #ifdef FASTQINDEXR_TIMING
    const auto plan_t0 = std::chrono::steady_clock::now();
    #endif
    const auto regions = buildDenseRegionsFromSortedUnique(
      local_ids,
      max_bridge_gap,
      max_region_records,
      bundle.record_size
      #ifdef FASTQINDEXR_TIMING
      , diag
      #endif
    );
    #ifdef FASTQINDEXR_TIMING
    if (diag != nullptr) {
      const auto plan_t1 = std::chrono::steady_clock::now();
      diag->time_region_plan_ms +=
        std::chrono::duration<double, std::milli>(plan_t1 - plan_t0).count();
    }
    #endif
    bool wrote_in_this_file = false;

    for (const auto& rg : regions) {
      const std::uint64_t region_start = rg.first;
      const std::uint64_t region_end = rg.second;
      const std::uint64_t line_start = region_start * static_cast<std::uint64_t>(bundle.record_size);
      const std::uint64_t line_count =
        (region_end - region_start + 1) * static_cast<std::uint64_t>(bundle.record_size);
      try {
        #ifdef FASTQINDEXR_TIMING
        if (diag != nullptr) {
          diag->extract_attempts++;
        }
        #endif
        const auto region_begin = std::lower_bound(local_ids.begin(), local_ids.end(), region_start);
        const auto region_end_it = std::upper_bound(local_ids.begin(), local_ids.end(), region_end);
        const std::uint64_t requested_in_region = static_cast<std::uint64_t>(region_end_it - region_begin);
        const std::uint64_t region_span = region_end - region_start + 1;
        const bool dense_region = (requested_in_region == region_span);
        if (dense_region) {
          #ifdef FASTQINDEXR_TIMING
          const auto dense_t0 = std::chrono::steady_clock::now();
          #endif
          std::uint64_t extracted_records = 0;
          extractor.extractRecords(
            Rcpp::as<std::string>(files[file_idx]),
            asGzipIndexedFile(*bundle.indexed_files[file_idx]).index.entries,
            line_start,
            line_count,
            bundle.record_size,
            region_start,
            bundle.format,
            need_qual,
            [&](std::uint64_t, const fastqindex_core::ExtractedRecord& rec_in) {
              ParsedRecord rec{rec_in.seq_id, rec_in.seq, rec_in.qual};
              writer->write(renderRecord(rec, resolved_output_type));
              wrote_in_this_file = true;
              n_written++;
              extracted_records++;
            }
          );
          if (extracted_records < region_span) {
            throw std::runtime_error("Partial extraction for requested region.");
          }
          if (extracted_records > region_span) {
            throw std::runtime_error("Malformed extracted line count for region.");
          }
          #ifdef FASTQINDEXR_TIMING
          if (diag != nullptr) {
            const auto dense_t1 = std::chrono::steady_clock::now();
            diag->time_indexed_dense_ms +=
              std::chrono::duration<double, std::milli>(dense_t1 - dense_t0).count();
          }
          #endif
        } else {
          #ifdef FASTQINDEXR_TIMING
          const auto sel_t0 = std::chrono::steady_clock::now();
          #endif
          const std::unordered_set<std::uint64_t> region_requested_set(region_begin, region_end_it);
          #ifdef FASTQINDEXR_TIMING
          fastqindex_core::SelectiveExtractTimings selective_timings;
          fastqindex_core::SelectiveExtractTimings* selective_timings_ptr =
            (diag != nullptr) ? &selective_timings : nullptr;
          #endif
          std::uint64_t selected_records = 0;
          extractor.extractSelectedRecords(
            Rcpp::as<std::string>(files[file_idx]),
            asGzipIndexedFile(*bundle.indexed_files[file_idx]).index.entries,
            line_start,
            line_count,
            bundle.record_size,
            region_requested_set,
            region_start,
            bundle.format,
            need_qual,
            [&](std::uint64_t, const fastqindex_core::ExtractedRecord& rec_in) {
              const ParsedRecord rec{rec_in.seq_id, rec_in.seq, rec_in.qual};
              writer->write(renderRecord(rec, resolved_output_type));
              wrote_in_this_file = true;
              n_written++;
              selected_records++;
            }
            #ifdef FASTQINDEXR_TIMING
            , selective_timings_ptr
            #endif
          );
          if (selected_records < requested_in_region) {
            throw std::runtime_error("Partial extraction for requested region.");
          }
          #ifdef FASTQINDEXR_TIMING
          if (diag != nullptr) {
            const auto sel_t1 = std::chrono::steady_clock::now();
            diag->time_indexed_selective_ms +=
              std::chrono::duration<double, std::milli>(sel_t1 - sel_t0).count();
            diag->time_selective_seek_init_ms += selective_timings_ptr->time_seek_init_ms;
            diag->time_selective_inflate_ms += selective_timings_ptr->time_inflate_ms;
            diag->time_selective_parse_ms += selective_timings_ptr->time_parse_ms;
            diag->time_selective_line_split_ms += selective_timings_ptr->time_line_split_ms;
            diag->time_selective_materialize_ms += selective_timings_ptr->time_materialize_ms;
            diag->time_selective_callback_ms += selective_timings_ptr->time_callback_ms;
          }
          #endif
        }
      } catch (const std::exception& e) {
        #ifdef FASTQINDEXR_TIMING
        recordExtractFailure(diag, e.what());
        #endif
        if (!wrote_in_this_file) {
          n_written += sequentialFallbackStreamToFile(
            Rcpp::as<std::string>(files[file_idx]),
            bundle.format,
            requested_set,
            need_qual,
            resolved_output_type,
            writer.get()
            #ifdef FASTQINDEXR_TIMING
            , diag
            #endif
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
            max_region_records,
            mode
            #ifdef FASTQINDEXR_TIMING
            , diag
            #endif
          );
        }
        break;
      } catch (...) {
        #ifdef FASTQINDEXR_TIMING
        recordExtractFailure(diag, "Unknown indexed extraction failure.");
        #endif
        if (!wrote_in_this_file) {
          n_written += sequentialFallbackStreamToFile(
            Rcpp::as<std::string>(files[file_idx]),
            bundle.format,
            requested_set,
            need_qual,
            resolved_output_type,
            writer.get()
            #ifdef FASTQINDEXR_TIMING
            , diag
            #endif
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
            max_region_records,
            mode
            #ifdef FASTQINDEXR_TIMING
            , diag
            #endif
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
  std::uint64_t max_region_records,
  ExtractExecutionMode mode
  #ifdef FASTQINDEXR_TIMING
  ,
  ExtractionDiagnostics* diag
  #endif
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
      max_region_records,
      mode
      #ifdef FASTQINDEXR_TIMING
      ,
      diag
      #endif
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
          false
          #ifdef FASTQINDEXR_TIMING
          ,
          diag
          #endif
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
      if (isPlainIndexedFile(*bundle.indexed_files[file_idx]) && mode == ExtractExecutionMode::Indexed) {
        const auto& offs = asPlainIndexedFile(*bundle.indexed_files[file_idx]).plain_header_byte_offsets;
        for (auto lid : local_ids) {
          if (lid >= offs.size()) {
            throw std::runtime_error("Requested id out of range for plain indexed file.");
          }
          ParsedRecord rec = readPlainRecordAtOffset(
            Rcpp::as<std::string>(files[file_idx]),
            offs[static_cast<size_t>(lid)],
            bundle.format,
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
        continue;
      }
      #ifdef FASTQINDEXR_TIMING
      const auto plan_t0 = std::chrono::steady_clock::now();
      #endif
      const auto regions = buildDenseRegionsFromSortedUnique(
        local_ids,
        max_bridge_gap,
        max_region_records,
        bundle.record_size
        #ifdef FASTQINDEXR_TIMING
        ,
        diag
        #endif
      );
      #ifdef FASTQINDEXR_TIMING
      if (diag != nullptr) {
        const auto plan_t1 = std::chrono::steady_clock::now();
        diag->time_region_plan_ms +=
          std::chrono::duration<double, std::milli>(plan_t1 - plan_t0).count();
      }
      #endif
      bool fell_back = false;

      for (const auto& rg : regions) {
        const std::uint64_t region_start = rg.first;
        const std::uint64_t region_end = rg.second;
        const std::uint64_t line_start = region_start * static_cast<std::uint64_t>(bundle.record_size);
        const std::uint64_t line_count =
          (region_end - region_start + 1) * static_cast<std::uint64_t>(bundle.record_size);
        try {
          #ifdef FASTQINDEXR_TIMING
          if (diag != nullptr) {
            diag->extract_attempts++;
          }
          #endif
          const auto region_begin = std::lower_bound(local_ids.begin(), local_ids.end(), region_start);
          const auto region_end_it = std::upper_bound(local_ids.begin(), local_ids.end(), region_end);
          const std::uint64_t requested_in_region = static_cast<std::uint64_t>(region_end_it - region_begin);
          const std::uint64_t region_span = region_end - region_start + 1;
          const bool dense_region = (requested_in_region == region_span);
          if (dense_region) {
            #ifdef FASTQINDEXR_TIMING
            const auto dense_t0 = std::chrono::steady_clock::now();
            #endif
            std::uint64_t extracted_records = 0;
            extractor.extractRecords(
              Rcpp::as<std::string>(files[file_idx]),
              asGzipIndexedFile(*bundle.indexed_files[file_idx]).index.entries,
              line_start,
              line_count,
              bundle.record_size,
              region_start,
              bundle.format,
              false,
              [&](std::uint64_t, const fastqindex_core::ExtractedRecord& rec_in) {
                const ParsedRecord rec{rec_in.seq_id, rec_in.seq, ""};
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
                extracted_records++;
              }
            );
            if (extracted_records < region_span) {
              throw std::runtime_error("Partial extraction for requested region.");
            }
            if (extracted_records > region_span) {
              throw std::runtime_error("Malformed extracted line count for region.");
            }
            #ifdef FASTQINDEXR_TIMING
            if (diag != nullptr) {
              const auto dense_t1 = std::chrono::steady_clock::now();
              diag->time_indexed_dense_ms +=
                std::chrono::duration<double, std::milli>(dense_t1 - dense_t0).count();
            }
            #endif
          } else {
            #ifdef FASTQINDEXR_TIMING
            const auto sel_t0 = std::chrono::steady_clock::now();
            #endif
            const std::unordered_set<std::uint64_t> region_requested_set(region_begin, region_end_it);
            #ifdef FASTQINDEXR_TIMING
            fastqindex_core::SelectiveExtractTimings selective_timings;
            fastqindex_core::SelectiveExtractTimings* selective_timings_ptr =
              (diag != nullptr) ? &selective_timings : nullptr;
            #endif
            std::uint64_t selected_records = 0;
            extractor.extractSelectedRecords(
              Rcpp::as<std::string>(files[file_idx]),
              asGzipIndexedFile(*bundle.indexed_files[file_idx]).index.entries,
              line_start,
              line_count,
              bundle.record_size,
              region_requested_set,
              region_start,
              bundle.format,
              false,
              [&](std::uint64_t, const fastqindex_core::ExtractedRecord& rec_in) {
                const ParsedRecord rec{rec_in.seq_id, rec_in.seq, ""};
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
                selected_records++;
              }
              #ifdef FASTQINDEXR_TIMING
              ,
              selective_timings_ptr
              #endif
            );
            if (selected_records < requested_in_region) {
              throw std::runtime_error("Partial extraction for requested region.");
            }
            #ifdef FASTQINDEXR_TIMING
            if (diag != nullptr) {
              const auto sel_t1 = std::chrono::steady_clock::now();
              diag->time_indexed_selective_ms +=
                std::chrono::duration<double, std::milli>(sel_t1 - sel_t0).count();
              diag->time_selective_seek_init_ms += selective_timings_ptr->time_seek_init_ms;
              diag->time_selective_inflate_ms += selective_timings_ptr->time_inflate_ms;
              diag->time_selective_parse_ms += selective_timings_ptr->time_parse_ms;
              diag->time_selective_line_split_ms += selective_timings_ptr->time_line_split_ms;
              diag->time_selective_materialize_ms += selective_timings_ptr->time_materialize_ms;
              diag->time_selective_callback_ms += selective_timings_ptr->time_callback_ms;
            }
            #endif
          }
        } catch (const std::exception& e) {
          #ifdef FASTQINDEXR_TIMING
          recordExtractFailure(diag, e.what());
          #endif
          fell_back = true;
          break;
        } catch (...) {
          #ifdef FASTQINDEXR_TIMING
          recordExtractFailure(diag, "Unknown indexed extraction failure.");
          #endif
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
          false
          #ifdef FASTQINDEXR_TIMING
          ,
          diag
          #endif
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

std::uint64_t countPlainTotalLines(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in.good()) {
    return 0;
  }
  std::uint64_t n = 0;
  std::string line;
  while (std::getline(in, line)) {
    n++;
  }
  return n;
}

SelectedRecordMap collectSequentialScanStandalone(
  Rcpp::CharacterVector files,
  const std::vector<std::uint64_t>& record_offsets,
  const std::string& format,
  const std::vector<std::uint64_t>& unique_sorted_global_ids,
  bool include_qual
  #ifdef FASTQINDEXR_TIMING
  ,
  ExtractionDiagnostics* diag
  #endif
) {
  const bool need_qual = (format == "fastq" && include_qual);
  SelectedRecordMap selected_global;
  selected_global.reserve(unique_sorted_global_ids.size());

  std::vector<std::vector<std::uint64_t>> local_ids_by_file(static_cast<size_t>(files.size()));
  for (auto gid : unique_sorted_global_ids) {
    auto up = std::upper_bound(record_offsets.begin(), record_offsets.end(), gid);
    if (up == record_offsets.begin() || up == record_offsets.end()) {
      throw std::runtime_error("Requested id out of range for scanned files.");
    }
    size_t file_idx = static_cast<size_t>((up - record_offsets.begin()) - 1);
    if (file_idx >= static_cast<size_t>(files.size())) {
      file_idx = static_cast<size_t>(files.size()) - 1;
    }
    const std::uint64_t local_id = gid - record_offsets[file_idx];
    local_ids_by_file[file_idx].push_back(local_id);
  }

  for (size_t file_idx = 0; file_idx < local_ids_by_file.size(); ++file_idx) {
    std::vector<std::uint64_t> local_ids = local_ids_by_file[file_idx];
    if (local_ids.empty()) {
      continue;
    }
    std::sort(local_ids.begin(), local_ids.end());
    local_ids.erase(std::unique(local_ids.begin(), local_ids.end()), local_ids.end());
    const std::unordered_set<std::uint64_t> requested_set(local_ids.begin(), local_ids.end());
    const std::string fpath = Rcpp::as<std::string>(files[static_cast<R_xlen_t>(file_idx)]);
    if (fileHasGzipMagic(fpath)) {
      scanSequentialRequestedRecords(
        fpath,
        format,
        need_qual,
        requested_set,
        [&](std::uint64_t local_id, const ParsedRecord& rec) {
          selected_global[record_offsets[file_idx] + local_id] = rec;
        }
        #ifdef FASTQINDEXR_TIMING
        ,
        diag
        #endif
      );
    } else {
      scanSequentialPlainRequestedRecords(
        fpath,
        format,
        need_qual,
        requested_set,
        [&](std::uint64_t local_id, const ParsedRecord& rec) {
          selected_global[record_offsets[file_idx] + local_id] = rec;
        }
        #ifdef FASTQINDEXR_TIMING
        ,
        diag
        #endif
      );
    }
  }

  return selected_global;
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
  const bool first_is_gzip = fileHasGzipMagic(paths.front());
  bundle.compression_kind = first_is_gzip ? CompressionKind::Gzip : CompressionKind::Plain;

  Rcpp::NumericVector per_file_counts(paths.size());
  Rcpp::CharacterVector file_compression(static_cast<R_xlen_t>(paths.size()));
  for (size_t i = 0; i < paths.size(); ++i) {
    const bool is_gzip = fileHasGzipMagic(paths[i]);
    if (is_gzip != first_is_gzip) {
      throw std::runtime_error("Mixed gzip/plain inputs are not supported in a single index.");
    }
    std::unique_ptr<IndexedFileBase> idxf;
    if (is_gzip) {
      std::unique_ptr<GzipIndexedFile> gzip_idx(new GzipIndexedFile());
      gzip_idx->index = indexer.createIndex(paths[i]);
      if (gzip_idx->index.total_lines % static_cast<std::uint64_t>(record_size) != 0) {
        throw std::runtime_error("Indexed line count is not a multiple of record size for: " + paths[i]);
      }
      gzip_idx->record_count = gzip_idx->index.total_lines / static_cast<std::uint64_t>(record_size);
      idxf = std::move(gzip_idx);
      file_compression[static_cast<R_xlen_t>(i)] = "gzip";
    } else {
      idxf = buildPlainFileIndex(paths[i], record_size);
      file_compression[static_cast<R_xlen_t>(i)] = "plain";
    }
    bundle.indexed_files.push_back(std::move(idxf));
    bundle.n_records += bundle.indexed_files.back()->record_count;
    bundle.record_offsets.push_back(bundle.n_records);
    per_file_counts[i] = static_cast<double>(bundle.indexed_files.back()->record_count);
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
    Rcpp::Named("file_compression") = file_compression,
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
  bundle.compression_kind = CompressionKind::Gzip;
  bundle.files = data_paths;
  bundle.record_offsets.push_back(0);

  Rcpp::NumericVector per_file_counts(data_paths.size());
  for (size_t i = 0; i < fqi_paths.size(); ++i) {
    std::unique_ptr<GzipIndexedFile> idxf = readSingleFqi(fqi_paths[i]);
    if (idxf->index.total_lines == 0 && std::filesystem::exists(data_paths[i])) {
      idxf->index.total_lines = countGzLines(data_paths[i]);
    }
    if (idxf->index.total_lines > 0 &&
        idxf->index.total_lines % static_cast<std::uint64_t>(record_size) != 0) {
      throw std::runtime_error(
        "Indexed line count is not a multiple of record size for: " + fqi_paths[i]
      );
    }
    idxf->record_count = idxf->index.total_lines / static_cast<std::uint64_t>(record_size);
    bundle.indexed_files.push_back(std::move(idxf));
    bundle.n_records += bundle.indexed_files.back()->record_count;
    bundle.record_offsets.push_back(bundle.n_records);
    per_file_counts[i] = static_cast<double>(bundle.indexed_files.back()->record_count);
  }

  Rcpp::XPtr<IndexBundle> index_ptr(new IndexBundle(std::move(bundle)), true);
  Rcpp::List payload = serializeIndexBundle(*index_ptr);

  Rcpp::NumericVector offsets(per_file_counts.size() + 1);
  for (R_xlen_t i = 0; i < offsets.size(); ++i) {
    offsets[i] = static_cast<double>((*index_ptr).record_offsets[static_cast<size_t>(i)]);
  }

  Rcpp::CharacterVector file_compression(static_cast<R_xlen_t>(data_paths.size()));
  for (R_xlen_t i = 0; i < file_compression.size(); ++i) {
    file_compression[i] = "gzip";
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
    Rcpp::Named("file_compression") = file_compression,
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
  double max_region_records,
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
  const std::uint64_t region_records = parsePositiveWhole(max_region_records, "max_region_records");
  const ExtractExecutionMode mode = parseExtractExecutionMode(extract_mode);
  #ifdef FASTQINDEXR_TIMING
  ExtractionDiagnostics diag;
  ExtractionDiagnostics* diag_ptr = diagnosticsPtrIfEnabled(&diag, diagnostics);
  #endif

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

  const bool need_parsed_qual = (bundle.format == "fastq" && include_qual);
  const SelectedRecordMap selected_global = collectRequestedRecords(
    files,
    bundle,
    unique_ids,
    need_parsed_qual,
    bridge_gap,
    region_records,
    mode
    #ifdef FASTQINDEXR_TIMING
    ,
    diag_ptr
    #endif
  );
  for (std::uint64_t gid : unique_ids) {
    auto rec_it = selected_global.find(gid);
    if (rec_it == selected_global.end()) {
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

  if (return_qual) {
    Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("seq_id") = out_id,
      Rcpp::Named("seq") = out_seq,
      Rcpp::Named("qual") = out_qual
    );
    #ifdef FASTQINDEXR_TIMING
    if (diagnostics && diag_ptr != nullptr) {
      out.attr("fastqindexr_diagnostics") = diagnosticsToList(diag);
    }
    #endif
    return out;
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("seq_id") = out_id,
    Rcpp::Named("seq") = out_seq
  );
  #ifdef FASTQINDEXR_TIMING
  if (diagnostics && diag_ptr != nullptr) {
    out.attr("fastqindexr_diagnostics") = diagnosticsToList(diag);
  }
  #endif
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
  double max_region_records,
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
  const std::uint64_t region_records = parsePositiveWhole(max_region_records, "max_region_records");
  const ExtractExecutionMode mode = parseExtractExecutionMode(extract_mode);
  #ifdef FASTQINDEXR_TIMING
  ExtractionDiagnostics diag;
  ExtractionDiagnostics* diag_ptr = diagnosticsPtrIfEnabled(&diag, diagnostics);
  #endif
  if (ids_zero_based.size() < 1) {
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("written") = 0.0);
    #ifdef FASTQINDEXR_TIMING
    if (diagnostics && diag_ptr != nullptr) {
      out.attr("fastqindexr_diagnostics") = diagnosticsToList(diag);
    }
    #endif
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
      region_records,
      mode
      #ifdef FASTQINDEXR_TIMING
      ,
      diag_ptr
      #endif
    );
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("written") = written);
    #ifdef FASTQINDEXR_TIMING
    if (diagnostics && diag_ptr != nullptr) {
      out.attr("fastqindexr_diagnostics") = diagnosticsToList(diag);
    }
    #endif
    return out;
  }

  const SelectedRecordMap selected_global = collectRequestedRecords(
    files,
    bundle,
    requested,
    need_parsed_qual,
    bridge_gap,
    region_records,
    mode
    #ifdef FASTQINDEXR_TIMING
    ,
    diag_ptr
    #endif
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
  #ifdef FASTQINDEXR_TIMING
  if (diagnostics && diag_ptr != nullptr) {
    out.attr("fastqindexr_diagnostics") = diagnosticsToList(diag);
  }
  #endif
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
  double max_region_records,
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
  const std::uint64_t region_records = parsePositiveWhole(max_region_records, "max_region_records");
  const ExtractExecutionMode mode = parseExtractExecutionMode(extract_mode);
  #ifdef FASTQINDEXR_TIMING
  ExtractionDiagnostics diag;
  ExtractionDiagnostics* diag_ptr = diagnosticsPtrIfEnabled(&diag, diagnostics);
  #endif
  if (renumber_mode != "none" &&
      renumber_mode != "zero_based" &&
      renumber_mode != "one_based") {
    throw std::runtime_error("Invalid renumber mode.");
  }
  if (ids_zero_based.size() == 0) {
    Rcpp::Environment biostrings = Rcpp::Environment::namespace_env("Biostrings");
    Rcpp::Function dna_string_set = biostrings["DNAStringSet"];
    SEXP empty = dna_string_set(Rcpp::CharacterVector(), Rcpp::Named("use.names") = true);
    #ifdef FASTQINDEXR_TIMING
    if (diagnostics && diag_ptr != nullptr) {
      Rf_setAttrib(empty, Rf_install("fastqindexr_diagnostics"), diagnosticsToList(diag));
    }
    #endif
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
    region_records,
    mode
    #ifdef FASTQINDEXR_TIMING
    ,
    diag_ptr
    #endif
  );
  #ifdef FASTQINDEXR_TIMING
  if (diagnostics && diag_ptr != nullptr) {
    Rf_setAttrib(out, Rf_install("fastqindexr_diagnostics"), diagnosticsToList(diag));
  }
  #endif
  return out;
}

// [[Rcpp::export]]
std::string cpp_detect_input_format(std::string path) {
  return openAndDetectType(path);
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_scan_record_offsets(Rcpp::CharacterVector files, std::string type) {
  if (type != "fasta" && type != "fastq") {
    throw std::runtime_error("`type` must be 'fasta' or 'fastq'.");
  }
  const int record_size = (type == "fastq") ? 4 : 2;
  Rcpp::NumericVector offsets(files.size() + 1);
  offsets[0] = 0.0;
  std::uint64_t cum = 0;
  for (R_xlen_t i = 0; i < files.size(); ++i) {
    const std::string p = Rcpp::as<std::string>(files[i]);
    std::uint64_t total_lines = 0;
    if (fileHasGzipMagic(p)) {
      total_lines = countGzLines(p);
    } else {
      total_lines = countPlainTotalLines(p);
    }
    if (total_lines % static_cast<std::uint64_t>(record_size) != 0) {
      throw std::runtime_error("Line count is not a multiple of record size for: " + p);
    }
    cum += total_lines / static_cast<std::uint64_t>(record_size);
    offsets[i + 1] = static_cast<double>(cum);
  }
  return offsets;
}

// [[Rcpp::export]]
Rcpp::List cpp_extract_sequences_streaming(
  Rcpp::CharacterVector files,
  std::string type,
  Rcpp::NumericVector ids_zero_based,
  Rcpp::NumericVector record_offsets_r,
  bool include_qual,
  bool diagnostics
) {
  if (type != "fasta" && type != "fastq") {
    throw std::runtime_error("`type` must be 'fasta' or 'fastq'.");
  }
  if (static_cast<size_t>(record_offsets_r.size()) != static_cast<size_t>(files.size()) + 1) {
    throw std::runtime_error("`record_offsets` must have length length(files) + 1.");
  }
  #ifdef FASTQINDEXR_TIMING
  ExtractionDiagnostics diag;
  ExtractionDiagnostics* diag_ptr = diagnosticsPtrIfEnabled(&diag, diagnostics);
  #endif

  const bool return_qual = (type == "fastq" && include_qual);
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

  std::vector<std::uint64_t> roff;
  roff.reserve(static_cast<size_t>(record_offsets_r.size()));
  for (R_xlen_t i = 0; i < record_offsets_r.size(); ++i) {
    roff.push_back(asUInt64(record_offsets_r[i], "record_offsets"));
  }

  const std::vector<std::uint64_t> requested = parseRequestedIds(ids_zero_based);
  R_xlen_t n = ids_zero_based.size();
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

  const SelectedRecordMap selected_global = collectSequentialScanStandalone(
    files,
    roff,
    type,
    unique_ids,
    include_qual
    #ifdef FASTQINDEXR_TIMING
    ,
    diag_ptr
    #endif
  );

  Rcpp::CharacterVector out_id(n);
  Rcpp::CharacterVector out_seq(n);
  Rcpp::CharacterVector out_qual;
  if (return_qual) {
    out_qual = Rcpp::CharacterVector(n);
  }

  for (std::uint64_t gid : unique_ids) {
    auto rec_it = selected_global.find(gid);
    if (rec_it == selected_global.end()) {
      throw std::runtime_error("Could not resolve all requested ids (streaming scan).");
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

  if (return_qual) {
    Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("seq_id") = out_id,
      Rcpp::Named("seq") = out_seq,
      Rcpp::Named("qual") = out_qual
    );
    #ifdef FASTQINDEXR_TIMING
    if (diagnostics && diag_ptr != nullptr) {
      out.attr("fastqindexr_diagnostics") = diagnosticsToList(diag);
    }
    #endif
    return out;
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("seq_id") = out_id,
    Rcpp::Named("seq") = out_seq
  );
  #ifdef FASTQINDEXR_TIMING
  if (diagnostics && diag_ptr != nullptr) {
    out.attr("fastqindexr_diagnostics") = diagnosticsToList(diag);
  }
  #endif
  return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_extract_sequences_to_file_streaming(
  Rcpp::CharacterVector files,
  std::string source_type,
  Rcpp::NumericVector ids_zero_based,
  Rcpp::NumericVector record_offsets_r,
  std::string output_type,
  std::string outfile,
  bool append,
  bool compress,
  bool diagnostics
) {
  if (source_type != "fasta" && source_type != "fastq") {
    throw std::runtime_error("`source_type` must be 'fasta' or 'fastq'.");
  }
  if (static_cast<size_t>(record_offsets_r.size()) != static_cast<size_t>(files.size()) + 1) {
    throw std::runtime_error("`record_offsets` must have length length(files) + 1.");
  }
  #ifdef FASTQINDEXR_TIMING
  ExtractionDiagnostics diag;
  ExtractionDiagnostics* diag_ptr = diagnosticsPtrIfEnabled(&diag, diagnostics);
  #endif

  if (ids_zero_based.size() < 1) {
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("written") = 0.0);
    #ifdef FASTQINDEXR_TIMING
    if (diagnostics && diag_ptr != nullptr) {
      out.attr("fastqindexr_diagnostics") = diagnosticsToList(diag);
    }
    #endif
    return out;
  }

  std::vector<std::uint64_t> roff;
  roff.reserve(static_cast<size_t>(record_offsets_r.size()));
  for (R_xlen_t i = 0; i < record_offsets_r.size(); ++i) {
    roff.push_back(asUInt64(record_offsets_r[i], "record_offsets"));
  }

  const std::string resolved_output_type = resolveOutputType(output_type, source_type);
  const std::vector<std::uint64_t> requested = parseRequestedIds(ids_zero_based);
  const bool need_parsed_qual = (resolved_output_type == "fastq");

  if (requestIsSortedUnique(requested) && !append) {
    std::unique_ptr<OutputWriter> writer = makeOutputWriter(outfile, false, compress);
    std::uint64_t n_written = 0;
    std::vector<std::vector<std::uint64_t>> local_ids_by_file(static_cast<size_t>(files.size()));
    for (auto gid : requested) {
      auto up = std::upper_bound(roff.begin(), roff.end(), gid);
      if (up == roff.begin() || up == roff.end()) {
        throw std::runtime_error("Requested id out of range for streaming extract.");
      }
      size_t file_idx = static_cast<size_t>((up - roff.begin()) - 1);
      if (file_idx >= static_cast<size_t>(files.size())) {
        file_idx = static_cast<size_t>(files.size()) - 1;
      }
      const std::uint64_t local_id = gid - roff[file_idx];
      local_ids_by_file[file_idx].push_back(local_id);
    }
    for (size_t file_idx = 0; file_idx < local_ids_by_file.size(); ++file_idx) {
      std::vector<std::uint64_t> local_ids = std::move(local_ids_by_file[file_idx]);
      if (local_ids.empty()) {
        continue;
      }
      const std::unordered_set<std::uint64_t> requested_set(local_ids.begin(), local_ids.end());
      const bool need_qual = (source_type == "fastq" && need_parsed_qual);
      const std::string fpath = Rcpp::as<std::string>(files[static_cast<R_xlen_t>(file_idx)]);
      if (fileHasGzipMagic(fpath)) {
        scanSequentialRequestedRecords(
          fpath,
          source_type,
          need_qual,
          requested_set,
          [&](std::uint64_t, const ParsedRecord& rec) {
            writer->write(renderRecord(rec, resolved_output_type));
            n_written++;
          }
          #ifdef FASTQINDEXR_TIMING
          ,
          diag_ptr
          #endif
        );
      } else {
        scanSequentialPlainRequestedRecords(
          fpath,
          source_type,
          need_qual,
          requested_set,
          [&](std::uint64_t, const ParsedRecord& rec) {
            writer->write(renderRecord(rec, resolved_output_type));
            n_written++;
          }
          #ifdef FASTQINDEXR_TIMING
          ,
          diag_ptr
          #endif
        );
      }
    }
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("written") = static_cast<double>(n_written));
    #ifdef FASTQINDEXR_TIMING
    if (diagnostics && diag_ptr != nullptr) {
      out.attr("fastqindexr_diagnostics") = diagnosticsToList(diag);
    }
    #endif
    return out;
  }

  std::vector<std::uint64_t> unique_sorted = requested;
  std::sort(unique_sorted.begin(), unique_sorted.end());
  unique_sorted.erase(std::unique(unique_sorted.begin(), unique_sorted.end()), unique_sorted.end());

  const SelectedRecordMap selected_global = collectSequentialScanStandalone(
    files,
    roff,
    source_type,
    unique_sorted,
    need_parsed_qual
    #ifdef FASTQINDEXR_TIMING
    ,
    diag_ptr
    #endif
  );

  std::unique_ptr<OutputWriter> writer = makeOutputWriter(outfile, append, compress);
  std::uint64_t written = 0;
  for (auto gid : requested) {
    auto it = selected_global.find(gid);
    if (it == selected_global.end()) {
      throw std::runtime_error("Could not resolve all requested ids (streaming extract).");
    }
    writer->write(renderRecord(it->second, resolved_output_type));
    written++;
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("written") = static_cast<double>(written)
  );
  #ifdef FASTQINDEXR_TIMING
  if (diagnostics && diag_ptr != nullptr) {
    out.attr("fastqindexr_diagnostics") = diagnosticsToList(diag);
  }
  #endif
  return out;
}

// [[Rcpp::export]]
SEXP cpp_extract_sequences_dnastringset_streaming(
  Rcpp::CharacterVector files,
  std::string source_type,
  Rcpp::NumericVector ids_zero_based,
  Rcpp::NumericVector record_offsets_r,
  double chunk_chars,
  bool diagnostics,
  std::string renumber_mode
) {
  if (source_type != "fasta" && source_type != "fastq") {
    throw std::runtime_error("`source_type` must be 'fasta' or 'fastq'.");
  }
  if (static_cast<size_t>(record_offsets_r.size()) != static_cast<size_t>(files.size()) + 1) {
    throw std::runtime_error("`record_offsets` must have length length(files) + 1.");
  }
  if (!std::isfinite(chunk_chars) || chunk_chars <= 0.0) {
    throw std::runtime_error("`chunk_chars` must be a finite positive value.");
  }
  if (renumber_mode != "none" &&
      renumber_mode != "zero_based" &&
      renumber_mode != "one_based") {
    throw std::runtime_error("Invalid renumber mode.");
  }
  #ifdef FASTQINDEXR_TIMING
  ExtractionDiagnostics diag;
  ExtractionDiagnostics* diag_ptr = diagnosticsPtrIfEnabled(&diag, diagnostics);
  #endif

  if (ids_zero_based.size() == 0) {
    Rcpp::Environment biostrings = Rcpp::Environment::namespace_env("Biostrings");
    Rcpp::Function dna_string_set = biostrings["DNAStringSet"];
    SEXP empty = dna_string_set(Rcpp::CharacterVector(), Rcpp::Named("use.names") = true);
    #ifdef FASTQINDEXR_TIMING
    if (diagnostics && diag_ptr != nullptr) {
      Rf_setAttrib(empty, Rf_install("fastqindexr_diagnostics"), diagnosticsToList(diag));
    }
    #endif
    return empty;
  }

  std::vector<std::uint64_t> roff;
  roff.reserve(static_cast<size_t>(record_offsets_r.size()));
  for (R_xlen_t i = 0; i < record_offsets_r.size(); ++i) {
    roff.push_back(asUInt64(record_offsets_r[i], "record_offsets"));
  }

  const std::vector<std::uint64_t> requested = parseRequestedIds(ids_zero_based);
  std::vector<std::uint64_t> unique_ids = requested;
  std::sort(unique_ids.begin(), unique_ids.end());
  unique_ids.erase(std::unique(unique_ids.begin(), unique_ids.end()), unique_ids.end());

  const SelectedRecordMap selected_global = collectSequentialScanStandalone(
    files,
    roff,
    source_type,
    unique_ids,
    false
    #ifdef FASTQINDEXR_TIMING
    ,
    diag_ptr
    #endif
  );

  std::vector<std::string> seq_chunk;
  std::vector<std::string> id_chunk;
  seq_chunk.reserve(1024);
  id_chunk.reserve(1024);
  Rcpp::List dna_chunks;
  std::uint64_t chunk_chars_u = static_cast<std::uint64_t>(chunk_chars);
  std::uint64_t chunk_chars_acc = 0;
  std::uint64_t output_index = 0;

  for (auto gid : requested) {
    auto it = selected_global.find(gid);
    if (it == selected_global.end()) {
      throw std::runtime_error("Could not resolve all requested ids (streaming DNA extract).");
    }
    appendRecordToDNAChunks(
      it->second,
      chunk_chars_u,
      &seq_chunk,
      &id_chunk,
      &chunk_chars_acc,
      &output_index,
      &dna_chunks,
      renumber_mode
    );
  }

  if (!seq_chunk.empty()) {
    const std::uint64_t chunk_start = output_index - seq_chunk.size();
    flushDNAChunk(&seq_chunk, &id_chunk, &dna_chunks, renumber_mode, chunk_start);
  }

  if (dna_chunks.size() == 0) {
    Rcpp::Environment biostrings = Rcpp::Environment::namespace_env("Biostrings");
    Rcpp::Function dna_string_set = biostrings["DNAStringSet"];
    SEXP empty = dna_string_set(Rcpp::CharacterVector(), Rcpp::Named("use.names") = true);
    #ifdef FASTQINDEXR_TIMING
    if (diagnostics && diag_ptr != nullptr) {
      Rf_setAttrib(empty, Rf_install("fastqindexr_diagnostics"), diagnosticsToList(diag));
    }
    #endif
    return empty;
  }
  if (dna_chunks.size() == 1) {
    SEXP one = dna_chunks[0];
    #ifdef FASTQINDEXR_TIMING
    if (diagnostics && diag_ptr != nullptr) {
      Rf_setAttrib(one, Rf_install("fastqindexr_diagnostics"), diagnosticsToList(diag));
    }
    #endif
    return one;
  }
  Rcpp::Environment base = Rcpp::Environment::base_env();
  Rcpp::Function do_call = base["do.call"];
  Rcpp::Function concat = base["c"];
  SEXP out = do_call(concat, dna_chunks);
  #ifdef FASTQINDEXR_TIMING
  if (diagnostics && diag_ptr != nullptr) {
    Rf_setAttrib(out, Rf_install("fastqindexr_diagnostics"), diagnosticsToList(diag));
  }
  #endif
  return out;
}
