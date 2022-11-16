#ifndef OPENMC_STATE_POINT_H
#define OPENMC_STATE_POINT_H

#include <cstdint>

#include "hdf5.h"

#include "openmc/capi.h"
#include "openmc/particle.h"
#include "openmc/vector.h"

#ifdef OPENMC_MCPL
#include <mcpl.h>
#endif

namespace openmc {

void load_state_point();
vector<int64_t> calculate_surf_source_size();
void write_source_point(const char* filename, bool surf_source_bank = false);
void write_source_bank(hid_t group_id, bool surf_source_bank);
void read_source_bank(
  hid_t group_id, vector<SourceSite>& sites, bool distribute);
void write_tally_results_nr(hid_t file_id);
void restart_set_keff();
void write_unstructured_mesh_results();

#ifdef OPENMC_MCPL
void write_mcpl_source_point(const char* filename, bool surf_source_bank = false);
void write_mcpl_source_bank(mcpl_outfile_t file_id, bool surf_source_bank);
#endif

} // namespace openmc
#endif // OPENMC_STATE_POINT_H
