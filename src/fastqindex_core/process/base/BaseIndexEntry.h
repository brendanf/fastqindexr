/**
 * Copyright (c) 2019 DKFZ - ODCF
 *
 * Distributed under the MIT License (license terms are at
 * https://github.com/dkfz-odcf/FastqIndEx/blob/master/LICENSE.txt).
 */

#ifndef FASTQINDEX_BASEINDEXENTRY_H
#define FASTQINDEX_BASEINDEXENTRY_H

namespace fastqindex_core {

struct BaseIndexEntry {
  bool operator==(const BaseIndexEntry&) const { return true; }
};

}  // namespace fastqindex_core

#endif //FASTQINDEX_BASEINDEXENTRY_H
