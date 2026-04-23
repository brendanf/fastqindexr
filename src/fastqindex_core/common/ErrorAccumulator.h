/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 */

#ifndef FASTQINDEX_ERRORACCUMULATOR_H
#define FASTQINDEX_ERRORACCUMULATOR_H

#include <sstream>
#include <string>
#include <vector>

namespace fastqindex_core {

class ErrorAccumulator {
  typedef const std::string& _cstr;

private:
  std::vector<std::string> errorMessages;
  static int verbosity;

public:
  virtual ~ErrorAccumulator() = default;

  static void setVerbosity(int verbosity);
  static bool verbosityIsSetToDebug();

  static void always(
    _cstr s0, _cstr s1 = "", _cstr s2 = "", _cstr s3 = "",
    _cstr s4 = "", _cstr s5 = "", _cstr s6 = ""
  );
  static void debug(
    _cstr s0, _cstr s1 = "", _cstr s2 = "", _cstr s3 = "",
    _cstr s4 = "", _cstr s5 = "", _cstr s6 = ""
  );
  static void info(const std::string& msg);
  static void warning(const std::string& msg);
  static void severe(const std::string& msg);

  virtual std::vector<std::string> getErrorMessages();

  void addErrorMessage(
    _cstr s0, _cstr s1 = "", _cstr s2 = "", _cstr s3 = "",
    _cstr s4 = "", _cstr s5 = "", _cstr s6 = ""
  );
  static std::string join(
    _cstr s0, _cstr s1 = "", _cstr s2 = "", _cstr s3 = "",
    _cstr s4 = "", _cstr s5 = "", _cstr s6 = ""
  );

  static std::vector<std::string> concatenateVectors(
    const std::vector<std::string>& l,
    const std::vector<std::string>& r
  );
  static std::vector<std::string> concatenateVectors(
    const std::vector<std::string>& a,
    const std::vector<std::string>& b,
    const std::vector<std::string>& c
  );
};

}  // namespace fastqindex_core

#endif //FASTQINDEX_ERRORACCUMULATOR_H
