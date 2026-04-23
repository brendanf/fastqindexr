/**
 * fastqindexr: R console sink for vendored ErrorAccumulator logging.
 *
 * Kept in this compilation unit so Rcpp stays out of the unity-built
 * fastqindex_core sources included from fastqindex_core_bridge.cpp.
 */

#include "fastqindexr_console.h"

#include <Rcpp.h>

void fastqindexr_log_to_r_console(const std::string& msg) {
  Rcpp::Rcerr << msg << "\n";
}
