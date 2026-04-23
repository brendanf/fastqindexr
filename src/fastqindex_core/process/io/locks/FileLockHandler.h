/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 */

#ifndef FASTQINDEX_FILELOCKHANDLER_H
#define FASTQINDEX_FILELOCKHANDLER_H

#include "process/io/locks/LockHandler.h"

#include <filesystem>
#include <mutex>
#include <stdio.h>

using namespace std;
using std::filesystem::path;

namespace fastqindex_core {

class FileLockHandler : public LockHandler {

private:
  mutex methodMutex;
  bool readLockActive = false;
  bool writeLockActive = false;
  FILE* lockedFileHandle = nullptr;

protected:
  path lockedFile;

public:
  explicit FileLockHandler(path file);
  ~FileLockHandler() override;
  path getLockedFile();
  bool readLock();
  bool writeLock();
  bool hasLock();
  void unlock();
};

}  // namespace fastqindex_core

#endif //FASTQINDEX_FILELOCKHANDLER_H
