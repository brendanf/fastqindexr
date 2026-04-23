/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 */

#include "FileSource.h"

namespace fastqindex_core {

FileSource::FileSource(const path& file)
  : file(IOHelper::fullPath(file)), lockHandler(file) {}

FileSource::~FileSource() {
  close();
}

bool FileSource::fulfillsPremises() {
  path check_path = file;
  bool is_valid = exists();
  if (is_valid) {
    if (isSymlink()) {
      check_path = read_symlink(check_path);
    }
    if (!is_regular_file(check_path)) {
      is_valid = false;
      addErrorMessage("'", toString(), "' does not point to a file.");
    }
  } else {
    addErrorMessage("File ", toString(), " does not exist.");
  }
  return is_valid;
}

bool FileSource::open() {
  if (!fStream.is_open()) {
    fStream.open(file, std::ifstream::binary);
  }
  return fStream.is_open();
}

bool FileSource::openWithReadLock() {
  if (!lockHandler.readLock()) {
    return false;
  }
  return open();
}

bool FileSource::close() {
  if (fStream.is_open()) {
    fStream.close();
  }
  if (lockHandler.hasLock()) {
    lockHandler.unlock();
  }
  return true;
}

int64_t FileSource::read(Bytef* targetBuffer, int numberOfBytes) {
  fStream.read(reinterpret_cast<char*>(targetBuffer), numberOfBytes);
  const int64_t amount_read = fStream.gcount();
  totalReadBytes += amount_read;
  return amount_read;
}

int FileSource::readChar() {
  Bytef result = 0;
  const int res = static_cast<int>(read(&result, 1));
  return res < 0 ? res : static_cast<int>(result);
}

int64_t FileSource::seek(int64_t nByte, bool absolute) {
  if (lastError()) {
    close();
    if (!open()) {
      return 0;
    }
  }

  if (absolute) {
    fStream.seekg(nByte, std::ifstream::beg);
  } else {
    fStream.seekg(nByte, std::ifstream::cur);
  }
  return (!fStream.fail() && !fStream.bad()) ? 1 : 0;
}

int64_t FileSource::skip(int64_t nBytes) {
  return seek(nBytes, false);
}

int64_t FileSource::tell() {
  if (fStream.is_open()) {
    return fStream.tellg();
  }
  return 0;
}

bool FileSource::canRead() {
  return tell() < size();
}

int FileSource::lastError() {
  return fStream.fail() || fStream.bad();
}

bool FileSource::isOpen() {
  return fStream.is_open();
}

bool FileSource::eof() {
  fStream.peek();
  return fStream.eof();
}

bool FileSource::isGood() {
  return fStream.good();
}

bool FileSource::empty() {
  return size() == 0;
}

bool FileSource::canWrite() {
  return false;
}

string FileSource::toString() {
  return file.string();
}

bool FileSource::hasLock() {
  return lockHandler.hasLock();
}

bool FileSource::unlock() {
  lockHandler.unlock();
  return !hasLock();
}

int64_t FileSource::rewind(int64_t nByte) {
  return seek(-nByte, false);
}

bool FileSource::isRegularFile() {
  if (isSymlink()) {
    return is_regular_file(read_symlink(file));
  }
  return is_regular_file(file);
}

}  // namespace fastqindex_core
