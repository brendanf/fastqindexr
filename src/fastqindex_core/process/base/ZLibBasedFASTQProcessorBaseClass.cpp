/**
 * Adapted from FastqIndEx/src/process/base/ZLibBasedFASTQProcessorBaseClass.cpp (MIT).
 *
 * fastqindexr note (omitted upstream functions):
 * - initializeZStreamForInflate() / initializeZStreamForRawInflate() wrappers.
 * - resetSlidingWindowIfNecessary(), decompressNextChunkOfData(), clearCurrentCompressedBlock()
 *   are implemented directly in Indexer/Extractor for simplified in-memory flow.
 */

#include "ZLibBasedFASTQProcessorBaseClass.h"

#include <cstring>

namespace fastqindex_core {

bool ZLibBasedFASTQProcessorBaseClass::initializeZStream(z_stream* strm, int mode) {
  std::memset(strm, 0, sizeof(z_stream));
  return inflateInit2(strm, mode) == Z_OK;
}

bool ZLibBasedFASTQProcessorBaseClass::readCompressedDataFromSource(
  std::istream& in,
  z_stream* strm,
  unsigned char* input,
  std::uint64_t* total_bytes_in
) {
  in.read(reinterpret_cast<char*>(input), CHUNK_SIZE);
  std::streamsize got = in.gcount();
  if (got <= 0) {
    return false;
  }
  strm->avail_in = static_cast<unsigned int>(got);
  strm->next_in = input;
  *total_bytes_in += static_cast<std::uint64_t>(got);
  return true;
}

bool ZLibBasedFASTQProcessorBaseClass::checkStreamForBlockEnd(const z_stream& strm) const {
  // fastqindexr change: same bit-test semantics as upstream helper.
  return (strm.data_type & 128) != 0 && !(strm.data_type & 64) != 0;
}

std::vector<std::string> ZLibBasedFASTQProcessorBaseClass::splitStr(const std::string& s) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == '\n') {
      out.push_back(cur);
      cur.clear();
    } else if (c != '\r') {
      cur.push_back(c);
    }
  }
  if (!cur.empty()) {
    out.push_back(cur);
  }
  return out;
}

}  // namespace fastqindex_core
