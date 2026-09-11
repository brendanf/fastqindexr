/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 *
 * fastqindexr (Windows):
 * - flock(2) / <sys/file.h> are POSIX and unavailable in Rtools. On
 *   Windows, holding an open FILE* is treated as the lock; this matches
 *   the single-process R usage (IndexReader via FileSource).
 * - fopen() needs a narrow path; path::c_str() is wchar_t* on Windows,
 *   so use path::string().
 */

#include "FileLockHandler.h"

#ifndef _WIN32
#include <sys/file.h>
#endif

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
  const string path_str = lockedFile.string();
  lockedFileHandle = fopen(path_str.c_str(), "rb");
  if (lockedFileHandle == nullptr) {
    return false;
  }
#ifdef _WIN32
  readLockActive = true;
  return true;
#else
  const bool result = flock(fileno(lockedFileHandle), LOCK_SH | LOCK_NB) == 0;
  if (result) {
    readLockActive = true;
  } else {
    fclose(lockedFileHandle);
    lockedFileHandle = nullptr;
  }
  return result;
#endif
}

bool FileLockHandler::writeLock() {
  lock_guard<mutex> lock(methodMutex);
  if (readLockActive || writeLockActive) {
    return false;
  }

  const string path_str = lockedFile.string();
  lockedFileHandle = fopen(path_str.c_str(), "wb");
  if (lockedFileHandle == nullptr) {
    return false;
  }

#ifdef _WIN32
  writeLockActive = true;
  return true;
#else
  const bool result = flock(fileno(lockedFileHandle), LOCK_EX | LOCK_NB) == 0;
  if (result) {
    writeLockActive = true;
  } else {
    fclose(lockedFileHandle);
    lockedFileHandle = nullptr;
  }
  return result;
#endif
}

bool FileLockHandler::hasLock() {
  return readLockActive || writeLockActive;
}

void FileLockHandler::unlock() {
  lock_guard<mutex> lock(methodMutex);
  readLockActive = false;
  writeLockActive = false;
  if (lockedFileHandle != nullptr) {
#ifndef _WIN32
    flock(fileno(lockedFileHandle), LOCK_UN);
#endif
    fclose(lockedFileHandle);
    lockedFileHandle = nullptr;
  }
}

}  // namespace fastqindex_core
