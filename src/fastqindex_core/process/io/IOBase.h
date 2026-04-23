/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 */

#ifndef FASTQINDEX_IOBASE_H
#define FASTQINDEX_IOBASE_H

#include "common/ErrorAccumulator.h"

#include <cstdint>
#include <string>

namespace fastqindex_core {

class IOBase : public ErrorAccumulator {

public:
  virtual bool fulfillsPremises() = 0;
  virtual bool open() = 0;
  virtual bool close() = 0;
  virtual bool hasLock() { return true; }
  virtual bool unlock() { return true; }
  virtual bool isOpen() = 0;
  virtual bool eof() = 0;
  virtual bool isGood() = 0;
  bool isBad() { return !isGood(); }
  virtual bool isFile() = 0;
  virtual bool isStream() = 0;
  virtual bool isSymlink() = 0;
  virtual bool exists() = 0;
  virtual int64_t size() = 0;
  virtual bool empty() = 0;
  virtual bool canRead() = 0;
  virtual bool canWrite() = 0;
  virtual int64_t seek(int64_t nByte, bool absolute) = 0;
  virtual int64_t skip(int64_t nByte) = 0;
  virtual int64_t rewind(int64_t nByte) { return seek(-nByte, false); }
  virtual std::string toString() = 0;
  virtual int64_t tell() = 0;
  virtual int lastError() = 0;
};

}  // namespace fastqindex_core

#endif //FASTQINDEX_IOBASE_H
