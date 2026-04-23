/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 */

#ifndef FASTQINDEX_INDEXREADER_H
#define FASTQINDEX_INDEXREADER_H

#include "common/CommonStructsAndConstants.h"
#include "common/ErrorAccumulator.h"
#include "process/base/IndexHeader.h"
#include "process/base/IndexEntry.h"
#include "process/io/Source.h"

#include <memory>
#include <vector>

namespace fastqindex_core {

class IndexReader : public ErrorAccumulator {

private:
  bool headerWasRead = false;
  bool readerIsOpen = false;
  std::shared_ptr<Source> indexFile;
  int64_t indicesLeft{0};
  int64_t indicesCount{0};
  IndexHeader readHeader;

  IndexHeader readIndexHeader();

public:
  explicit IndexReader(const std::shared_ptr<Source>& indexFile);
  ~IndexReader() override;

  bool tryOpenAndReadHeader();
  std::vector<std::shared_ptr<IndexEntry>> readIndexFile();
  std::vector<std::shared_ptr<IndexEntryV1>> readIndexFileV1();
  std::shared_ptr<IndexEntry> readIndexEntry();
  std::shared_ptr<IndexEntryV1> readIndexEntryV1();

  IndexHeader getIndexHeader() { return readHeader; }
  int64_t getIndicesLeft() { return indicesLeft; }
};

}  // namespace fastqindex_core

#endif //FASTQINDEX_INDEXREADER_H
