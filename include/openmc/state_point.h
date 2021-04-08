#ifndef OPENMC_STATE_POINT_H
#define OPENMC_STATE_POINT_H

#include <cstdint>
#include <vector>

#include "hdf5.h"

#include "openmc/capi.h"
#include "openmc/particle.h"

namespace openmc {

void load_state_point();
std::vector<int64_t> calculate_surf_source_size();
void write_source_point(const char* filename, bool surf_source_bank = false);
void write_source_bank(hid_t group_id, bool surf_source_bank);
void read_source_bank(hid_t group_id, std::vector<Particle::Bank>& sites, bool distribute);
void write_tally_results_nr(hid_t file_id);
void restart_set_keff();
void write_unstructured_mesh_results();

} // namespace openmc
#endif // OPENMC_STATE_POINT_H
