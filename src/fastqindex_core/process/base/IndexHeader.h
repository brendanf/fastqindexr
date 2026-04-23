/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 */

#ifndef FASTQINDEX_INDEXHEADER_H
#define FASTQINDEX_INDEXHEADER_H

#include "common/CommonStructsAndConstants.h"
#include "IndexEntry.h"
#include "IndexEntryV1.h"

namespace fastqindex_core {

struct IndexHeader {
  std::uint32_t indexWriterVersion{0};
  std::uint32_t sizeOfIndexEntry{0};
  std::uint32_t magicNumber = MAGIC_NUMBER;
  std::uint32_t blockInterval{0};
  int64_t numberOfEntries{0};
  int64_t linesInIndexedFile{0};
  bool dictionariesAreCompressed{false};
  Bytef placeholder[7]{0};
  int64_t reserved[59]{0};

  explicit IndexHeader(
    std::uint32_t binaryVersion,
    std::uint32_t sizeOfIndexEntry,
    std::uint32_t blockInterval,
    bool dictionariesAreCompressed
  ) {
    this->indexWriterVersion = binaryVersion;
    this->sizeOfIndexEntry = sizeOfIndexEntry;
    this->blockInterval = blockInterval;
    this->dictionariesAreCompressed = dictionariesAreCompressed;
  }

  IndexHeader() = default;

  bool operator==(const IndexHeader& rhs) const {
    return indexWriterVersion == rhs.indexWriterVersion &&
      sizeOfIndexEntry == rhs.sizeOfIndexEntry &&
      magicNumber == rhs.magicNumber &&
      blockInterval == rhs.blockInterval &&
      numberOfEntries == rhs.numberOfEntries &&
      linesInIndexedFile == rhs.linesInIndexedFile &&
      dictionariesAreCompressed == rhs.dictionariesAreCompressed;
  }

  bool operator!=(const IndexHeader& rhs) const { return !(rhs == *this); }
  explicit operator bool() const { return magicNumber == MAGIC_NUMBER; }
  bool operator!() { return !this->operator bool(); }
};

}  // namespace fastqindex_core

#endif //FASTQINDEX_INDEXHEADER_H
