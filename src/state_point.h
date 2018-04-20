#ifndef STATE_POINT_H
#define STATE_POINT_H

#include <cstdint>

#include "hdf5.h"
#include "openmc.h"

namespace openmc {

extern "C" void write_source_bank(hid_t group_id, int64_t* work_index,
                                  const struct Bank* bank);

} // namespace openmc
#endif // STATE_POINT_H
