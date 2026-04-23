/**
 * fastqindexr: bridge from vendored C++ to the R console.
 *
 * Declares a single logging hook implemented in fastqindexr_console.cpp so
 * vendored units (e.g. ErrorAccumulator.cpp pulled into
 * fastqindex_core_bridge.cpp) can avoid including Rcpp.
 */

#ifndef FASTQINDEXR_CONSOLE_H
#define FASTQINDEXR_CONSOLE_H

#include <string>

void fastqindexr_log_to_r_console(const std::string& msg);

#endif
