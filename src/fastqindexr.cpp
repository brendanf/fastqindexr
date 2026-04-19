#include <Rcpp.h>
#include <zlib.h>

#include "fastqindex_core/process/extract/Extractor.h"
#include "fastqindex_core/process/index/Indexer.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

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

std::unordered_map<std::string, IndexBundle> g_index_registry;
std::atomic<std::uint64_t> g_token_counter{1};

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

std::vector<std::pair<std::uint64_t, std::uint64_t>> denseRegions(std::vector<std::uint64_t> ids) {
  std::vector<std::pair<std::uint64_t, std::uint64_t>> regions;
  if (ids.empty()) {
    return regions;
  }
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
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

ParsedRecord parseRecordFromFixedLines(const std::string& type, const std::vector<std::string>& lines, size_t offset) {
  if (type == "fastq") {
    if (offset + 3 >= lines.size()) {
      throw std::runtime_error("Malformed FASTQ lines in extracted region.");
    }
    return ParsedRecord{trimPrefix(lines[offset], '@'), lines[offset + 1], lines[offset + 3]};
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
  std::unordered_map<std::uint64_t, ParsedRecord>* selected_global
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
        (*selected_global)[global_offset + local_id] = ParsedRecord{trimPrefix(a, '@'), b, d};
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

  const std::string token = "idx_" + std::to_string(g_token_counter.fetch_add(1));
  g_index_registry[token] = std::move(bundle);

  Rcpp::NumericVector offsets(per_file_counts.size() + 1);
  for (R_xlen_t i = 0; i < offsets.size(); ++i) {
    offsets[i] = static_cast<double>(g_index_registry[token].record_offsets[i]);
  }

  return Rcpp::List::create(
    Rcpp::Named("format") = type,
    Rcpp::Named("files") = files,
    Rcpp::Named("file_info") = fileInfoDataFrame(files),
    Rcpp::Named("file_record_counts") = per_file_counts,
    Rcpp::Named("file_record_offsets") = offsets,
    Rcpp::Named("n_records") = static_cast<double>(g_index_registry[token].n_records),
    Rcpp::Named("index_token") = token,
    Rcpp::Named("record_size") = record_size,
    Rcpp::Named("engine") = "upstream_adapted_v1"
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_extract_sequences(
  Rcpp::CharacterVector files,
  std::string type,
  Rcpp::NumericVector ids_zero_based,
  std::string index_token
) {
  auto it_bundle = g_index_registry.find(index_token);
  if (it_bundle == g_index_registry.end()) {
    throw std::runtime_error("Unknown/expired index token; recreate index.");
  }
  IndexBundle& bundle = it_bundle->second;
  if (bundle.format != type) {
    throw std::runtime_error("Type mismatch between index and request.");
  }
  if (static_cast<size_t>(files.size()) != bundle.files.size()) {
    throw std::runtime_error("Override `file` must contain same file count as index.");
  }

  if (ids_zero_based.size() < 1) {
    return Rcpp::List::create(
      Rcpp::Named("seq_id") = Rcpp::CharacterVector(),
      Rcpp::Named("seq") = Rcpp::CharacterVector(),
      Rcpp::Named("qual") = Rcpp::CharacterVector()
    );
  }

  std::vector<std::uint64_t> requested(ids_zero_based.size());
  for (R_xlen_t i = 0; i < ids_zero_based.size(); ++i) {
    requested[i] = static_cast<std::uint64_t>(ids_zero_based[i]);
  }

  std::unordered_map<std::uint64_t, ParsedRecord> selected_global;
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
    auto local_ids = local_ids_by_file[file_idx];
    if (local_ids.empty()) {
      continue;
    }
    std::sort(local_ids.begin(), local_ids.end());
    local_ids.erase(std::unique(local_ids.begin(), local_ids.end()), local_ids.end());
    std::unordered_set<std::uint64_t> requested_set(local_ids.begin(), local_ids.end());
    auto regions = denseRegions(local_ids);
    for (const auto& rg : regions) {
      std::uint64_t region_start = rg.first;
      std::uint64_t region_end = rg.second;
      std::uint64_t line_start = region_start * static_cast<std::uint64_t>(bundle.record_size);
      std::uint64_t line_count = (region_end - region_start + 1) * static_cast<std::uint64_t>(bundle.record_size);

      try {
        std::vector<std::string> lines = extractor.extract(
          Rcpp::as<std::string>(files[file_idx]),
          bundle.indexed_files[file_idx].index.entries,
          line_start,
          line_count
        );
        std::size_t extracted_records = lines.size() / static_cast<std::size_t>(bundle.record_size);
        for (std::size_t r = 0; r < extracted_records; ++r) {
          std::uint64_t local_id = region_start + static_cast<std::uint64_t>(r);
          if (requested_set.find(local_id) == requested_set.end()) {
            continue;
          }
          ParsedRecord rec = parseRecordFromFixedLines(
            bundle.format,
            lines,
            r * static_cast<std::size_t>(bundle.record_size)
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
          &selected_global
        );
      }
    }
  }

  R_xlen_t n = ids_zero_based.size();
  Rcpp::CharacterVector out_id(n);
  Rcpp::CharacterVector out_seq(n);
  Rcpp::CharacterVector out_qual(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    std::uint64_t gid = requested[i];
    auto it = selected_global.find(gid);
    if (it == selected_global.end()) {
      throw std::runtime_error("Could not resolve all requested ids with provided index/files.");
    }
    out_id[i] = it->second.id;
    out_seq[i] = it->second.seq;
    out_qual[i] = it->second.qual;
  }

  return Rcpp::List::create(
    Rcpp::Named("seq_id") = out_id,
    Rcpp::Named("seq") = out_seq,
    Rcpp::Named("qual") = out_qual
  );
}
