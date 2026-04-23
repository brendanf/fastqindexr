/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 */

#ifndef FASTQINDEX_FILESOURCE_H
#define FASTQINDEX_FILESOURCE_H

#include "common/IOHelper.h"
#include "process/io/Source.h"
#include "process/io/locks/FileLockHandler.h"

#include <filesystem>
#include <fstream>
#include <memory>

namespace fastqindex_core {

class FileSource : public Source {
private:
  path file;
  std::ifstream fStream;
  FileLockHandler lockHandler;

public:
  static shared_ptr<Source> from(const path& file) {
    return make_shared<FileSource>(file);
  }

  bool hasLock() override;
  bool unlock() override;
  int64_t rewind(int64_t nByte) override;

  explicit FileSource(const path& file);
  ~FileSource() override;

  bool fulfillsPremises() override;
  bool open() override;
  bool openWithReadLock() override;
  bool close() override;
  bool exists() override { return std::filesystem::exists(file); }
  bool isSymlink() override { return is_symlink(symlink_status(file)); }
  bool isRegularFile();
  int64_t size() override { return exists() ? file_size(file) : 0; }
  string absolutePath() { return file.string(); }
  bool isFile() override { return true; }
  bool isStream() override { return false; }
  int64_t read(Bytef* targetBuffer, int numberOfBytes) override;
  int readChar() override;
  int64_t seek(int64_t nByte, bool absolute) override;
  int64_t skip(int64_t nBytes) override;
  int64_t tell() override;
  bool canRead() override;
  int lastError() override;
  path getPath() { return file; }
  bool isOpen() override;
  bool eof() override;
  bool isGood() override;
  bool empty() override;
  bool canWrite() override;
  string toString() override;
};

}  // namespace fastqindex_core

#include "Source.h"

#endif //FASTQINDEX_FILESOURCE_H
