/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 */

#ifndef FASTQINDEX_SOURCE_H
#define FASTQINDEX_SOURCE_H

#include "process/io/IOBase.h"

#include <filesystem>
#include <memory>

using namespace std;
using namespace std::filesystem;

namespace fastqindex_core {

class Source : public IOBase {
protected:
  int64_t totalReadBytes{0};
  int64_t readStart{0};

  Source() = default;

public:
  ~Source() override = default;

  virtual bool openWithReadLock() = 0;
  virtual int64_t getTotalReadBytes() { return totalReadBytes; }

  virtual void setReadStart(int64_t startBytes) {
    readStart = startBytes;
  }

  virtual int64_t read(Bytef* targetBuffer, int numberOfBytes) = 0;
  virtual int readChar() = 0;
};

}  // namespace fastqindex_core

#endif //FASTQINDEX_SOURCE_H
