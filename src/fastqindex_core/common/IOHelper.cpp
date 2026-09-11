/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 *
 * fastqindexr (R / R CMD check):
 * - IOHelper::report: no std::cerr fallback when errorAccumulator is null
 *   (use a non-null accumulator; null is unsupported, message dropped; see
 *   IOHelper.h and ErrorAccumulator + fastqindexr_console.cpp for console
 *   policy).
 * - getApplicationPath: readlink(2) return value is checked, buffer is
 *   null-terminated, and failed reads return an empty path. This avoids
 *   GCC -Wunused-result, which R CMD check surfaces as a significant
 *   install-time WARNING on typical Linux toolchains.
 * - Windows: upstream uses POSIX pwd/unistd APIs that Rtools/MinGW does
 *   not provide. Guarded replacements keep the Unix bodies unchanged.
 */

#include "IOHelper.h"

#include <cerrno>
#include <cstdlib>
#include <fstream>
#ifdef _WIN32
#include <direct.h>
#include <fcntl.h>
#include <io.h>
#include <sys/stat.h>
#else
#include <pwd.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace fastqindex_core {

recursive_mutex IOHelper::iohelper_mtx;

void IOHelper::report(const stringstream& sstream, ErrorAccumulator* errorAccumulator) {
  if (errorAccumulator != nullptr) {
    errorAccumulator->addErrorMessage(sstream.str());
  }
  // fastqindexr: null accumulator — intentionally silent; see file header.
}

path IOHelper::getUserHomeDirectory() {
#ifdef _WIN32
  // <pwd.h> / getpwuid() are POSIX. USERPROFILE is the usual home dir;
  // HOME is a common MinGW/MSYS fallback.
  const char *home = std::getenv("USERPROFILE");
  if (home == nullptr || home[0] == '\0')
  {
    home = std::getenv("HOME");
  }
  if (home == nullptr || home[0] == '\0')
  {
    return path();
  }
  return path(string(home));
#else
  struct passwd* pw = getpwuid(getuid());
  if (pw == nullptr) {
    return path();
  }
  return path(string(pw->pw_dir));
#endif
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
#ifdef _WIN32
  // access()/W_OK are POSIX; CRT _access() uses mode 2 for write.
  auto not_writable = [](const path &p)
  {
    return _access(p.string().c_str(), 2) != 0;
  };
#else
  auto not_writable = [](const path &p)
  {
    return access(p.string().c_str(), W_OK) != 0;
  };
#endif
  if (exists(file) && not_writable(file))
  {
    sstream << "The '" << fileType << "' file '" << file.string()
            << "' exists but is not writeable.";
    error = true;
  }
  else if (exists(file.parent_path()) && not_writable(file.parent_path()))
  {
    sstream << "The parent folder for '" << fileType << "' file '" << file.string()
            << "' exists but is not writeable.";
    error = true;
  }
  else if (!exists(file.parent_path()))
  {
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
#ifdef _WIN32
  // mkdtemp(3) is POSIX. CRT _mktemp() fills the trailing X's; _mkdir
  // then creates that directory.
  path tmpl = temp_directory_path() / (prefix + "_XXXXXXXXXXXXXX");
  const string test_dir = tmpl.string();
  unique_ptr<char[]> buf(new char[test_dir.size() + 1]{0});
  test_dir.copy(buf.get(), test_dir.size(), 0);
  if (_mktemp(buf.get()) == nullptr)
  {
    return {false, path()};
  }
  if (_mkdir(buf.get()) != 0)
  {
    return {false, path()};
  }
  return {true, path(string(buf.get()))};
#else
  const auto temp_dir = temp_directory_path();
  const string test_dir = temp_dir.string() + "/" + prefix + "_XXXXXXXXXXXXXX";
  unique_ptr<char[]> buf(new char[test_dir.size() + 1]{0});
  test_dir.copy(buf.get(), test_dir.size(), 0);
  char* result = mkdtemp(buf.get());
  return {result != nullptr, path(result == nullptr ? "" : string(result))};
#endif
}

tuple<bool, path> IOHelper::createTempFile(const string& prefix) {
  lock_guard<recursive_mutex> lock_guard(iohelper_mtx);
#ifdef _WIN32
  // mkstemp(3) is POSIX. _mktemp() + exclusive _open() is the CRT
  // counterpart.
  path tmpl = temp_directory_path() / (prefix + "_XXXXXXXXXXXXXX");
  const string tmp_file = tmpl.string();
  unique_ptr<char[]> buf(new char[tmp_file.size() + 1]{0});
  tmp_file.copy(buf.get(), tmp_file.size(), 0);
  if (_mktemp(buf.get()) == nullptr)
  {
    return {false, path()};
  }
  const int result = _open(
      buf.get(), _O_CREAT | _O_EXCL | _O_RDWR, _S_IREAD | _S_IWRITE);
  if (result != -1)
  {
    _close(result);
  }
  return {result != -1, path(string(buf.get()))};
#else
  const auto temp_dir = temp_directory_path();
  const string tmp_file = temp_dir.string() + "/" + prefix + "_XXXXXXXXXXXXXX";
  unique_ptr<char[]> buf(new char[tmp_file.size() + 1]{0});
  tmp_file.copy(buf.get(), tmp_file.size(), 0);
  const int result = mkstemp(buf.get());
  if (result != -1) {
    close(result);
  }
  return {result != -1, path(string(buf.get()))};
#endif
}

tuple<bool, path> IOHelper::createTempFifo(const string& prefix) {
  lock_guard<recursive_mutex> lock_guard(iohelper_mtx);
#ifdef _WIN32
  // mkfifo(3) has no Windows equivalent in this form. Unused in the R
  // path; fail instead of creating a regular file.
  (void)prefix;
  return {false, path()};
#else
  auto [success, fifo_path] = createTempFile(prefix);
  remove(fifo_path);
  if (success) {
    mkfifo(fifo_path.string().c_str(), 0600);
  }
  return {success, fifo_path};
#endif
}

path IOHelper::fullPath(const path& file) {
  char buf[32768]{0};
#ifdef _WIN32
  // realpath(3) is POSIX; _fullpath() is the Windows CRT equivalent.
  if (_fullpath(buf, file.string().c_str(), sizeof(buf)) == nullptr)
  {
    return file;
  }
#else
  if (realpath(file.string().c_str(), buf) == nullptr) {
    return file;
  }
#endif
  return path(string(buf));
}

// readlink is annotated warn_unused_result on glibc; ignoring it made
// R CMD check fail the install step. We cap the read with sizeof(buf)-1,
// write an explicit NUL after n bytes, and return {} when readlink fails.
path IOHelper::getApplicationPath() {
#ifdef _WIN32
  // readlink("/proc/self/exe") is Linux-specific. Unused in the R path.
  return path();
#else
  char buf[32768]{0};
  ssize_t n = -1;
  if (exists("/proc")) {
    n = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
  } else if (exists("/user")) {
    n = readlink("/user/self/exe", buf, sizeof(buf) - 1);
  }
  if (n <= 0) {
    return path();
  }
  buf[static_cast<size_t>(n)] = '\0';
  return path(string(buf));
#endif
}

shared_ptr<unordered_map<string, string>> IOHelper::loadIniFile(
  const path&, const string&
) {
  return make_shared<unordered_map<string, string>>();
}

}  // namespace fastqindex_core
