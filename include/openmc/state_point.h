#ifndef OPENMC_STATE_POINT_H
#define OPENMC_STATE_POINT_H

#include <cstdint>

#include "hdf5.h"

#include "openmc/capi.h"

namespace openmc {

void load_state_point();
void write_source_point(const char* filename);
void write_source_bank(hid_t group_id);
void read_source_bank(hid_t group_id);
void write_tally_results_nr(hid_t file_id);
void restart_set_keff();

#ifdef DAGMC
void write_unstructured_mesh_results();
#endif

} // namespace openmc
#endif // OPENMC_STATE_POINT_H
