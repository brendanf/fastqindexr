/**
 * Adapted from FastqIndEx/src/common/CommonStructsAndConstants.h (MIT).
 *
 * fastqindexr note (omitted upstream items that lived in the upstream .cpp):
 * - FQI_BINARY / S3HELPER_BINARY globals are not included (no CLI helpers).
 * - MAGIC_NUMBER and storage-unit constants are not included (no .fqi binary IO).
 * - stoui() is not included (no CLI parsing path).
 */

#ifndef FASTQINDEXR_COMMONSTRUCTSANDCONSTANTS_H
#define FASTQINDEXR_COMMONSTRUCTSANDCONSTANTS_H

namespace fastqindex_core {

constexpr unsigned int WINDOW_SIZE = 32768U;
constexpr unsigned int CHUNK_SIZE = 16384U;

}  // namespace fastqindex_core

#endif
