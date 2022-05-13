#include "openmc/track_output.h"

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/position.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/vector.h"

#include "xtensor/xtensor.hpp"
#include <fmt/core.h>

#include <cstddef> // for size_t
#include <string>

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

//==============================================================================
// Non-member functions
//==============================================================================

void add_particle_track(Particle& p)
{
  p.tracks().emplace_back();
  p.tracks().back().particle = p.type();
}

void write_particle_track(Particle& p)
{
  p.tracks().back().states.push_back(p.get_track_state());
}

hid_t trackstate_type()
{
  // Create compound type for Position
  hid_t postype = H5Tcreate(H5T_COMPOUND, sizeof(struct Position));
  H5Tinsert(postype, "x", HOFFSET(Position, x), H5T_NATIVE_DOUBLE);
  H5Tinsert(postype, "y", HOFFSET(Position, y), H5T_NATIVE_DOUBLE);
  H5Tinsert(postype, "z", HOFFSET(Position, z), H5T_NATIVE_DOUBLE);

  // Create compound type for TrackState
  hid_t tracktype = H5Tcreate(H5T_COMPOUND, sizeof(struct TrackState));
  H5Tinsert(tracktype, "r", HOFFSET(TrackState, r), postype);
  H5Tinsert(tracktype, "u", HOFFSET(TrackState, u), postype);
  H5Tinsert(tracktype, "E", HOFFSET(TrackState, E), H5T_NATIVE_DOUBLE);
  H5Tinsert(tracktype, "time", HOFFSET(TrackState, time), H5T_NATIVE_DOUBLE);
  H5Tinsert(tracktype, "wgt", HOFFSET(TrackState, wgt), H5T_NATIVE_DOUBLE);
  H5Tinsert(tracktype, "cell_id", HOFFSET(TrackState, cell_id), H5T_NATIVE_INT);
  H5Tinsert(
    tracktype, "material_id", HOFFSET(TrackState, material_id), H5T_NATIVE_INT);

  H5Tclose(postype);
  return tracktype;
}

void finalize_particle_track(Particle& p)
{
  std::string filename =
    fmt::format("{}track_{}_{}_{}.h5", settings::path_output,
      simulation::current_batch, simulation::current_gen, p.id());

  // Determine number of coordinates for each particle
  vector<int> offsets;
  vector<int> particles;
  vector<TrackState> tracks;
  int offset = 0;
  for (auto& track_i : p.tracks()) {
    offsets.push_back(offset);
    particles.push_back(static_cast<int>(track_i.particle));
    offset += track_i.states.size();
    tracks.insert(tracks.end(), track_i.states.begin(), track_i.states.end());
  }
  offsets.push_back(offset);

#pragma omp critical(FinalizeParticleTrack)
  {
    hid_t file_id = file_open(filename, 'w');
    write_attribute(file_id, "filetype", "track");
    write_attribute(file_id, "version", VERSION_TRACK);
    write_attribute(file_id, "n_particles", p.tracks().size());
    write_attribute(file_id, "offsets", offsets);
    write_attribute(file_id, "particles", particles);

    // Create HDF5 datatype for TrackState
    hid_t track_type = trackstate_type();

    // Write array of TrackState to file
    hsize_t dims[] {static_cast<hsize_t>(tracks.size())};
    hid_t dspace = H5Screate_simple(1, dims, nullptr);
    hid_t dset = H5Dcreate(file_id, "tracks", track_type, dspace, H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, track_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, tracks.data());

    // Free resources
    H5Dclose(dset);
    H5Sclose(dspace);
    H5Tclose(track_type);
    file_close(file_id);
  }

  // Clear particle tracks
  p.tracks().clear();
}

} // namespace openmc
