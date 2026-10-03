#ifndef OPENMC_STATE_POINT_H
#define OPENMC_STATE_POINT_H

#include <cstdint>
#include <string>

#include "hdf5.h"

#include "openmc/capi.h"
#include "openmc/particle.h"
#include "openmc/shared_array.h"
#include "openmc/span.h"
#include "openmc/vector.h"

namespace openmc {

void load_state_point();

// By passing in a filename, source bank, and list of source indices
// on each MPI rank, this writes an HDF5 file which contains that
// information which can later be read in by read_source_bank
// (defined below). If you're writing code to write out a new kind
// of particle bank, this function is the one you want to use!
//
// For example, this is used to write both the surface source sites
// or fission source sites for eigenvalue continuation runs.
//
// This function ends up calling write_source_bank, and is responsible
// for opening the file to be written to and controlling whether the
// write is done in parallel (if compiled with parallel HDF5).
//
// bank_index is an exclusive parallel scan of the source_bank.size()
// values on each rank, used to create global indexing. This vector
// can be created by calling calculate_parallel_index_vector on
// source_bank.size() if such a vector is not already available.
//
// The source_bank variable is used as work space if MPI is used,
// so it cannot be given as a const span.
//
// If birth_mesh_tags is true, the fission bank birth tags of the sites are
// written as a "birth_mesh_bin" dataset together with a "birth_mesh_id"
// attribute. This should only be requested for a bank of fission sites banked
// in an eigenvalue calculation with a <birth_mesh> (see settings::birth_mesh).
void write_h5_source_point(const char* filename, span<SourceSite> source_bank,
  const vector<int64_t>& bank_index, bool birth_mesh_tags = false);

void write_source_point(std::string, span<SourceSite> source_bank,
  const vector<int64_t>& bank_index, bool use_mcpl,
  bool birth_mesh_tags = false);

// This appends a source bank specification to an HDF5 file
// that's already open. It is used internally by write_source_point. If
// birth_mesh_tags is true, the "birth_mesh_bin" dataset is also written.
void write_source_bank(hid_t group_id, span<SourceSite> source_bank,
  const vector<int64_t>& bank_index, bool birth_mesh_tags = false);

void read_source_bank(
  hid_t group_id, vector<SourceSite>& sites, bool distribute);
void write_tally_results_nr(hid_t file_id);
void restart_set_keff();
void write_unstructured_mesh_results();

} // namespace openmc
#endif // OPENMC_STATE_POINT_H
