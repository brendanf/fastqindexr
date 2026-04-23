/**
 * Adapted from FastqIndEx/src/common/CommonStructsAndConstants.h (MIT).
 *
 * fastqindexr note:
 * - Keeps upstream-like constants in one place for both in-memory indexing and
 *   vendored `.fqi` index reading.
 * - Omits CLI binary path globals and parsing helpers not needed in-package.
 */

#ifndef FASTQINDEXR_COMMONSTRUCTSANDCONSTANTS_H
#define FASTQINDEXR_COMMONSTRUCTSANDCONSTANTS_H

#include <cstdint>
#include <zlib.h>

namespace fastqindex_core {

constexpr unsigned int WINDOW_SIZE = 32768U;
constexpr unsigned int CHUNK_SIZE = 16384U;
constexpr unsigned int CLEAN_WINDOW_SIZE = WINDOW_SIZE + 1U;
constexpr std::uint32_t MAGIC_NUMBER = 0x04030201U;
constexpr std::int64_t kB = 1024;
constexpr std::int64_t MB = kB * 1024;
constexpr std::int64_t GB = MB * 1024;
constexpr std::int64_t TB = GB * 1024;
constexpr int DEFAULT_RECORD_SIZE = 4;

}  // namespace fastqindex_core

#endif
