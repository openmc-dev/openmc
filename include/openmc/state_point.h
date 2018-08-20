#ifndef OPENMC_STATE_POINT_H
#define OPENMC_STATE_POINT_H

#include <cstdint>

#include "hdf5.h"

#include "openmc/capi.h"

namespace openmc {

extern "C" void write_source_bank(hid_t group_id, int64_t* work_index,
                                  Bank* source_bank);
extern "C" void read_source_bank(hid_t group_id, int64_t* work_index,
                                 Bank* source_bank);

} // namespace openmc
#endif // OPENMC_STATE_POINT_H
