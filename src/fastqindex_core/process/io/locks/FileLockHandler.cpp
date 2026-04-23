/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 */

#include "FileLockHandler.h"

#include <sys/file.h>

namespace fastqindex_core {

FileLockHandler::FileLockHandler(path file) : lockedFile(std::move(file)) {}

FileLockHandler::~FileLockHandler() {
  unlock();
}

path FileLockHandler::getLockedFile() {
  return lockedFile;
}

bool FileLockHandler::readLock() {
  lock_guard<mutex> lock(methodMutex);
  if (readLockActive) {
    return true;
  }
  lockedFileHandle = fopen(lockedFile.c_str(), "rb");
  if (lockedFileHandle == nullptr) {
    return false;
  }
  const bool result = flock(fileno(lockedFileHandle), LOCK_SH | LOCK_NB) == 0;
  if (result) {
    readLockActive = true;
  } else {
    fclose(lockedFileHandle);
    lockedFileHandle = nullptr;
  }
  return result;
}

bool FileLockHandler::writeLock() {
  lock_guard<mutex> lock(methodMutex);
  if (readLockActive || writeLockActive) {
    return false;
  }

  lockedFileHandle = fopen(lockedFile.c_str(), "wb");
  if (lockedFileHandle == nullptr) {
    return false;
  }

  const bool result = flock(fileno(lockedFileHandle), LOCK_EX | LOCK_NB) == 0;
  if (result) {
    writeLockActive = true;
  } else {
    fclose(lockedFileHandle);
    lockedFileHandle = nullptr;
  }
  return result;
}

bool FileLockHandler::hasLock() {
  return readLockActive || writeLockActive;
}

void FileLockHandler::unlock() {
  lock_guard<mutex> lock(methodMutex);
  readLockActive = false;
  writeLockActive = false;
  if (lockedFileHandle != nullptr) {
    flock(fileno(lockedFileHandle), LOCK_UN);
    fclose(lockedFileHandle);
    lockedFileHandle = nullptr;
  }
}

}  // namespace fastqindex_core
