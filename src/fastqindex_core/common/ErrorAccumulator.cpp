/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 *
 * fastqindexr: stderr logging replaced with fastqindexr_log_to_r_console()
 * (fastqindexr_console.cpp) so this vendored TU stays free of Rcpp while
 * satisfying R's compiled-code rules for console output.
 */

#include "ErrorAccumulator.h"

#include "fastqindexr_console.h"

namespace fastqindex_core {

int ErrorAccumulator::verbosity = 0;

void ErrorAccumulator::setVerbosity(int level) {
  if (level >= 0 && level <= 3) {
    verbosity = level;
  }
}

bool ErrorAccumulator::verbosityIsSetToDebug() {
  return verbosity >= 3;
}

void ErrorAccumulator::always(
  _cstr s0, _cstr s1, _cstr s2, _cstr s3, _cstr s4, _cstr s5, _cstr s6
) {
  fastqindexr_log_to_r_console(ErrorAccumulator::join(s0, s1, s2, s3, s4, s5, s6));
}

void ErrorAccumulator::debug(
  _cstr s0, _cstr s1, _cstr s2, _cstr s3, _cstr s4, _cstr s5, _cstr s6
) {
  if (verbosityIsSetToDebug()) {
    fastqindexr_log_to_r_console(
      ErrorAccumulator::join(s0, s1, s2, s3, s4, s5, s6)
    );
  }
}

void ErrorAccumulator::info(const std::string& msg) {
  if (verbosity >= 2) {
    fastqindexr_log_to_r_console(msg);
  }
}

void ErrorAccumulator::warning(const std::string& msg) {
  if (verbosity >= 1) {
    fastqindexr_log_to_r_console(msg);
  }
}

void ErrorAccumulator::severe(const std::string& msg) {
  fastqindexr_log_to_r_console(msg);
}

std::vector<std::string> ErrorAccumulator::getErrorMessages() {
  return errorMessages;
}

void ErrorAccumulator::addErrorMessage(
  _cstr s0, _cstr s1, _cstr s2, _cstr s3, _cstr s4, _cstr s5, _cstr s6
) {
  errorMessages.emplace_back(join(s0, s1, s2, s3, s4, s5, s6));
}

std::string ErrorAccumulator::join(
  _cstr s0, _cstr s1, _cstr s2, _cstr s3, _cstr s4, _cstr s5, _cstr s6
) {
  std::ostringstream stream;
  stream << s0 << s1 << s2 << s3 << s4 << s5 << s6;
  return stream.str();
}

std::vector<std::string> ErrorAccumulator::concatenateVectors(
  const std::vector<std::string>& l,
  const std::vector<std::string>& r
) {
  std::vector<std::string> merged;
  merged.reserve(l.size() + r.size());
  merged.insert(merged.end(), l.begin(), l.end());
  merged.insert(merged.end(), r.begin(), r.end());
  return merged;
}

std::vector<std::string> ErrorAccumulator::concatenateVectors(
  const std::vector<std::string>& a,
  const std::vector<std::string>& b,
  const std::vector<std::string>& c
) {
  std::vector<std::string> merged;
  merged.reserve(a.size() + b.size() + c.size());
  merged.insert(merged.end(), a.begin(), a.end());
  merged.insert(merged.end(), b.begin(), b.end());
  merged.insert(merged.end(), c.begin(), c.end());
  return merged;
}

}  // namespace fastqindex_core
