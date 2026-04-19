/**
 * Adapted from FastqIndEx/src/process/base/ZLibBasedFASTQProcessorBaseClass.h (MIT).
 */

#ifndef FASTQINDEXR_ZLIBBASEDFASTQPROCESSORBASECLASS_H
#define FASTQINDEXR_ZLIBBASEDFASTQPROCESSORBASECLASS_H

#include "../../common/CommonStructsAndConstants.h"

#include <cstdint>
#include <istream>
#include <string>
#include <vector>

#include <zlib.h>

namespace fastqindex_core {

class ZLibBasedFASTQProcessorBaseClass {
 protected:
  bool initializeZStream(z_stream* strm, int mode);
  bool readCompressedDataFromSource(
    std::istream& in,
    z_stream* strm,
    unsigned char* input,
    std::uint64_t* total_bytes_in
  );
  bool checkStreamForBlockEnd(const z_stream& strm) const;

 public:
  static std::vector<std::string> splitStr(const std::string& s);
};

}  // namespace fastqindex_core

#endif
