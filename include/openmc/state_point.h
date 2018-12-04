#ifndef OPENMC_STATE_POINT_H
#define OPENMC_STATE_POINT_H

#include <cstdint>

#include "hdf5.h"

#include "openmc/capi.h"

namespace openmc {

void write_source_point(const char* filename);
extern "C" void write_source_bank(hid_t group_id);
extern "C" void read_source_bank(hid_t group_id);
extern "C" void write_tally_results_nr(hid_t file_id);
extern "C" void restart_set_keff();

} // namespace openmc
#endif // OPENMC_STATE_POINT_H
