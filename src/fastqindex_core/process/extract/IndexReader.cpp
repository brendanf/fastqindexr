/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 */

#include "IndexReader.h"

#include <cstring>

namespace {
constexpr unsigned int kIndexWriterVersion = 1U;
}

namespace fastqindex_core {

IndexReader::IndexReader(const std::shared_ptr<Source>& indexFile) {
  this->indexFile = indexFile;
}

IndexReader::~IndexReader() {
  if (this->indexFile) {
    this->indexFile->close();
  }
}

bool IndexReader::tryOpenAndReadHeader() {
  if (readerIsOpen) {
    return true;
  }

  if (!indexFile->exists()) {
    addErrorMessage("Index file '", indexFile->toString(), "' does not exist.");
    return false;
  }
  const bool gotLock = indexFile->openWithReadLock();
  if (!gotLock) {
    addErrorMessage("Could not get a lock for index file '", indexFile->toString(), "'.");
    return false;
  }

  const int64_t fileSize = static_cast<int64_t>(indexFile->size());
  const int64_t headerSize = static_cast<int64_t>(sizeof(IndexHeader));
  if (headerSize > fileSize) {
    addErrorMessage(
      "Index file '", indexFile->toString(),
      "' is smaller than the minimum size and cannot be read."
    );
    indexFile->close();
    return false;
  }

  if (this->indexFile->isBad()) {
    addErrorMessage("Could not open input stream for index file '", indexFile->toString(), "'.");
    indexFile->close();
    return false;
  }

  this->readIndexHeader();

  unsigned int sizeOfIndexEntry = 0;
  if (this->readHeader.indexWriterVersion == 1) {
    sizeOfIndexEntry = sizeof(IndexEntryV1);
  } else {
    addErrorMessage(
      "Index version '",
      std::to_string(this->readHeader.indexWriterVersion),
      "' is not readable (max supported ",
      std::to_string(kIndexWriterVersion),
      ")."
    );
    indexFile->close();
    return false;
  }

  if (this->readHeader.dictionariesAreCompressed) {
    if (fileSize <= headerSize) {
      addErrorMessage("Compressed index file does not contain entries.");
      indexFile->close();
      return false;
    }
  } else if (headerSize + static_cast<int64_t>(sizeOfIndexEntry) > fileSize) {
    addErrorMessage("Index file is smaller than minimum non-compressed size.");
    indexFile->close();
    return false;
  }

  if (!this->readHeader.dictionariesAreCompressed &&
      0 != (fileSize - headerSize) % static_cast<int64_t>(sizeOfIndexEntry)) {
    addErrorMessage("Stored index version and content size mismatch.");
    indexFile->close();
    return false;
  }

  if (this->readHeader.numberOfEntries > 0) {
    this->indicesLeft = this->readHeader.numberOfEntries;
  } else if (!this->readHeader.dictionariesAreCompressed) {
    this->indicesLeft = (fileSize - headerSize) / static_cast<int64_t>(sizeOfIndexEntry);
  }

  if (this->indicesLeft == 0) {
    addErrorMessage("Could not determine entry count in index file '", indexFile->toString(), "'.");
    indexFile->close();
    return false;
  }

  this->indicesCount = indicesLeft;
  this->readerIsOpen = true;
  return readerIsOpen;
}

IndexHeader IndexReader::readIndexHeader() {
  IndexHeader header;
  indexFile->read(reinterpret_cast<Bytef*>(&header), sizeof(IndexHeader));
  this->readHeader = header;
  headerWasRead = true;
  return header;
}

std::vector<std::shared_ptr<IndexEntry>> IndexReader::readIndexFile() {
  std::vector<std::shared_ptr<IndexEntry>> convertedLines;
  if (!tryOpenAndReadHeader()) {
    addErrorMessage("Could not read index file '", indexFile->toString(), "' during open.");
    return convertedLines;
  }

  if (this->readHeader.indexWriterVersion == 1) {
    const auto entries = readIndexFileV1();
    for (const auto& entry : entries) {
      convertedLines.emplace_back(entry->toIndexEntry());
    }
  }
  return convertedLines;
}

std::vector<std::shared_ptr<IndexEntryV1>> IndexReader::readIndexFileV1() {
  std::vector<std::shared_ptr<IndexEntryV1>> convertedLines;
  while (indicesLeft > 0) {
    auto indexEntry = readIndexEntryV1();
    if (indexEntry) {
      convertedLines.emplace_back(indexEntry);
    }
  }
  return convertedLines;
}

std::shared_ptr<IndexEntry> IndexReader::readIndexEntry() {
  if (!tryOpenAndReadHeader()) {
    addErrorMessage("Could not read index file '", indexFile->toString(), "' during open.");
    return std::shared_ptr<IndexEntry>(nullptr);
  }

  if (this->readHeader.indexWriterVersion == 1) {
    return readIndexEntryV1()->toIndexEntry();
  }

  addErrorMessage("BUG: Reached impossible code path in IndexReader::readIndexEntry()");
  return std::shared_ptr<IndexEntry>(nullptr);
}

std::shared_ptr<IndexEntryV1> IndexReader::readIndexEntryV1() {
  if (!readerIsOpen) {
    addErrorMessage("BUG: open IndexReader via tryOpenAndReadHeader() first.");
    return std::shared_ptr<IndexEntryV1>(nullptr);
  }

  if (indexFile->eof() || indicesLeft <= 0) {
    addErrorMessage("Stream finished and no entries are left.");
    return std::shared_ptr<IndexEntryV1>(nullptr);
  }

  auto entry = std::make_shared<IndexEntryV1>();
  const int headerSize = sizeof(IndexEntryV1) - sizeof(entry->dictionary);
  indexFile->read(reinterpret_cast<Bytef*>(entry.get()), headerSize);
  if (entry->compressedDictionarySize == 0) {
    indexFile->read(reinterpret_cast<Bytef*>(entry.get()) + headerSize, sizeof(entry->dictionary));
  } else {
    std::memset(reinterpret_cast<void*>(entry->dictionary), 0, WINDOW_SIZE);
    indexFile->read(
      reinterpret_cast<Bytef*>(entry.get()) + headerSize,
      entry->compressedDictionarySize
    );
  }
  indicesLeft--;
  return entry;
}

}  // namespace fastqindex_core
