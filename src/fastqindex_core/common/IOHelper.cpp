/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 */

#include "IOHelper.h"

#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <pwd.h>
#include <sys/stat.h>
#include <unistd.h>

namespace fastqindex_core {

recursive_mutex IOHelper::iohelper_mtx;

void IOHelper::report(const stringstream& sstream, ErrorAccumulator* errorAccumulator) {
  if (errorAccumulator != nullptr) {
    errorAccumulator->addErrorMessage(sstream.str());
  } else {
    cerr << sstream.str() << "\n";
  }
}

path IOHelper::getUserHomeDirectory() {
  struct passwd* pw = getpwuid(getuid());
  if (pw == nullptr) {
    return path();
  }
  return path(string(pw->pw_dir));
}

bool IOHelper::checkFileReadability(
  const path& file, const string& fileType, ErrorAccumulator* errorAccumulator
) {
  bool is_valid = exists(file);
  if (!is_valid) {
    stringstream sstream;
    sstream << "The " << fileType << " file '" << file.string()
            << "' could not be found or is inaccessible.";
    report(sstream, errorAccumulator);
    return false;
  }

  ifstream file_stream(file);
  is_valid = file_stream.good();
  file_stream.close();
  if (!is_valid) {
    stringstream sstream;
    sstream << "The " << fileType << " file '" << file.string()
            << "' could not be read.";
    report(sstream, errorAccumulator);
  }
  return is_valid;
}

bool IOHelper::checkFileWriteability(
  const path& file, const string& fileType, ErrorAccumulator* errorAccumulator
) {
  stringstream sstream;
  bool error = false;
  if (exists(file) && access(file.string().c_str(), W_OK) != 0) {
    sstream << "The '" << fileType << "' file '" << file.string()
            << "' exists but is not writeable.";
    error = true;
  } else if (exists(file.parent_path()) &&
             access(file.parent_path().string().c_str(), W_OK) != 0) {
    sstream << "The parent folder for '" << fileType << "' file '" << file.string()
            << "' exists but is not writeable.";
    error = true;
  } else if (!exists(file.parent_path())) {
    sstream << "The parent folder for '" << fileType << "' file '" << file.string()
            << "' does not exist.";
    error = true;
  }

  if (error) {
    report(sstream, errorAccumulator);
  }
  return !error;
}

tuple<bool, path> IOHelper::createTempDir(const string& prefix) {
  lock_guard<recursive_mutex> lock_guard(iohelper_mtx);
  const auto temp_dir = temp_directory_path();
  const string test_dir = temp_dir.string() + "/" + prefix + "_XXXXXXXXXXXXXX";
  unique_ptr<char[]> buf(new char[test_dir.size() + 1]{0});
  test_dir.copy(buf.get(), test_dir.size(), 0);
  char* result = mkdtemp(buf.get());
  return {result != nullptr, path(result == nullptr ? "" : string(result))};
}

tuple<bool, path> IOHelper::createTempFile(const string& prefix) {
  lock_guard<recursive_mutex> lock_guard(iohelper_mtx);
  const auto temp_dir = temp_directory_path();
  const string tmp_file = temp_dir.string() + "/" + prefix + "_XXXXXXXXXXXXXX";
  unique_ptr<char[]> buf(new char[tmp_file.size() + 1]{0});
  tmp_file.copy(buf.get(), tmp_file.size(), 0);
  const int result = mkstemp(buf.get());
  if (result != -1) {
    close(result);
  }
  return {result != -1, path(string(buf.get()))};
}

tuple<bool, path> IOHelper::createTempFifo(const string& prefix) {
  lock_guard<recursive_mutex> lock_guard(iohelper_mtx);
  auto [success, fifo_path] = createTempFile(prefix);
  remove(fifo_path);
  if (success) {
    mkfifo(fifo_path.string().c_str(), 0600);
  }
  return {success, fifo_path};
}

path IOHelper::fullPath(const path& file) {
  char buf[32768]{0};
  if (realpath(file.string().c_str(), buf) == nullptr) {
    return file;
  }
  return path(string(buf));
}

path IOHelper::getApplicationPath() {
  char buf[32768]{0};
  if (exists("/proc")) {
    readlink("/proc/self/exe", buf, 32768);
  } else if (exists("/user")) {
    readlink("/user/self/exe", buf, 32768);
  }
  return path(string(buf));
}

shared_ptr<unordered_map<string, string>> IOHelper::loadIniFile(
  const path&, const string&
) {
  return make_shared<unordered_map<string, string>>();
}

}  // namespace fastqindex_core
