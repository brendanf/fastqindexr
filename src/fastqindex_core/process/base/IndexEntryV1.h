/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 */

#ifndef FASTQINDEX_INDEXENTRYV1_H
#define FASTQINDEX_INDEXENTRYV1_H

#include "common/CommonStructsAndConstants.h"
#include "process/base/BaseIndexEntry.h"
#include "IndexEntry.h"

#include <cstring>
#include <memory>
#include <stdexcept>
#include <zlib.h>

namespace fastqindex_core {

struct IndexEntryV1;
typedef std::shared_ptr<IndexEntryV1> IndexEntryV1_S;

struct IndexEntryV1 : public BaseIndexEntry {
  std::uint64_t blockIndex{0};
  std::uint64_t blockOffsetInRawFile{0};
  std::uint64_t startingLineInEntry{0};
  std::uint32_t offsetToNextLineStart{0};
  unsigned char bits{0};
  unsigned char reserved{0};
  std::uint16_t compressedDictionarySize{0};
  Bytef dictionary[WINDOW_SIZE]{0};

  static IndexEntryV1_S from(
    unsigned char bits,
    std::uint64_t blockID,
    std::uint32_t offsetOfFirstValidLine,
    std::uint64_t offsetInRawFile,
    std::uint64_t startingLineInEntry
  ) {
    return std::make_shared<IndexEntryV1>(
      bits, blockID, offsetOfFirstValidLine, offsetInRawFile, startingLineInEntry
    );
  }

  IndexEntryV1(
    unsigned char bits,
    std::uint64_t blockID,
    std::uint32_t offsetOfFirstValidLine,
    std::uint64_t offsetInRawFile,
    std::uint64_t startingLineInEntry
  ) :
    blockIndex(blockID),
    blockOffsetInRawFile(offsetInRawFile),
    startingLineInEntry(startingLineInEntry),
    offsetToNextLineStart(offsetOfFirstValidLine),
    bits(bits) {}

  IndexEntryV1() = default;

  bool operator==(const IndexEntryV1& rhs) const {
    return bits == rhs.bits &&
      blockIndex == rhs.blockIndex &&
      offsetToNextLineStart == rhs.offsetToNextLineStart &&
      blockOffsetInRawFile == rhs.blockOffsetInRawFile &&
      startingLineInEntry == rhs.startingLineInEntry;
  }

  bool operator!=(const IndexEntryV1& rhs) const { return !(rhs == *this); }

  std::shared_ptr<IndexEntry> toIndexEntry() {
    auto indexLine = std::make_shared<IndexEntry>();
    indexLine->block_index = blockIndex;
    indexLine->block_offset_in_raw_file = blockOffsetInRawFile;
    indexLine->starting_line_in_entry = startingLineInEntry;
    indexLine->offset_to_next_line_start = offsetToNextLineStart;
    indexLine->bits = bits;
    indexLine->dictionary.resize(WINDOW_SIZE);
    if (compressedDictionarySize > 0) {
      z_stream strm;
      std::memset(&strm, 0, sizeof(strm));
      strm.next_in = const_cast<Bytef*>(dictionary);
      strm.avail_in = compressedDictionarySize;
      strm.next_out = indexLine->dictionary.data();
      strm.avail_out = WINDOW_SIZE;
      if (inflateInit(&strm) != Z_OK) {
        throw std::runtime_error("inflateInit failed for compressed dictionary.");
      }
      const int ret = inflate(&strm, Z_FINISH);
      inflateEnd(&strm);
      if (ret != Z_STREAM_END || strm.total_out != WINDOW_SIZE) {
        throw std::runtime_error("Failed to inflate compressed dictionary.");
      }
    } else {
      std::memcpy(indexLine->dictionary.data(), dictionary, WINDOW_SIZE);
    }
    return indexLine;
  }
};

}  // namespace fastqindex_core

#endif //FASTQINDEX_INDEXENTRYV1_H
