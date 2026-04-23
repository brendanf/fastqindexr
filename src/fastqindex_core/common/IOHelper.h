/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 */

#ifndef FASTQINDEX_IOHELPER_H
#define FASTQINDEX_IOHELPER_H

#include "ErrorAccumulator.h"

#include <filesystem>
#include <memory>
#include <mutex>
#include <string>
#include <tuple>
#include <unordered_map>

using namespace std;
using namespace std::filesystem;

namespace fastqindex_core {

/**
 * fastqindexr: For report() and the check* helpers, pass a non-null
 * ErrorAccumulator; a null pointer drops messages (no stderr fallback).
 */
class IOHelper {

private:
  static recursive_mutex iohelper_mtx;

public:
  static void report(const stringstream& sstream, ErrorAccumulator* errorAccumulator);
  static path getUserHomeDirectory();
  static tuple<bool, path> createTempDir(const string& prefix);
  static tuple<bool, path> createTempFile(const string& prefix);
  static tuple<bool, path> createTempFifo(const string& prefix);
  static bool checkFileReadability(const path& file, const string& fileType, ErrorAccumulator* errorAccumulator);
  static bool checkFileWriteability(const path& file, const string& fileType, ErrorAccumulator* errorAccumulator);
  static path fullPath(const path& file);
  static path getApplicationPath();

  // Upstream has an INI parser dependency; it is not needed for index reading.
  static shared_ptr<unordered_map<string, string>> loadIniFile(
    const path& /* file */, const string& /* section */
  );
};

}  // namespace fastqindex_core

#endif //FASTQINDEX_IOHELPER_H
