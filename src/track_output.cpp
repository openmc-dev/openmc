#include "openmc/track_output.h"

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/position.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"

#include "xtensor/xtensor.hpp"

#include <cstddef> // for size_t
#include <sstream>
#include <string>
#include <vector>

// Explicit vector template specialization of threadprivate variable outside of
// the openmc namespace for the picky Intel compiler.
template class std::vector<std::vector<openmc::Position>>;

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

// Forward declaration needed in order to declare tracks as threadprivate
extern std::vector<std::vector<Position>> tracks;
#pragma omp threadprivate(tracks)

std::vector<std::vector<Position>> tracks;

//==============================================================================
// Non-member functions
//==============================================================================

void add_particle_track()
{
  tracks.emplace_back();
}

void write_particle_track(const Particle& p)
{
  tracks.back().push_back(p.r());
}

void finalize_particle_track(const Particle& p)
{
  std::stringstream filename;
  filename << settings::path_output << "track_" << simulation::current_batch
    << '_' << simulation::current_gen << '_' << p.id_ << ".h5";

  // Determine number of coordinates for each particle
  std::vector<int> n_coords;
  for (auto& coords : tracks) {
    n_coords.push_back(coords.size());
  }

  #pragma omp critical (FinalizeParticleTrack)
  {
    hid_t file_id = file_open(filename.str().c_str(), 'w');
    write_attribute(file_id, "filetype", "track");
    write_attribute(file_id, "version", VERSION_TRACK);
    write_attribute(file_id, "n_particles", tracks.size());
    write_attribute(file_id, "n_coords", n_coords);
    for (int i = 1; i <= tracks.size(); ++i) {
      const auto& t {tracks[i-1]};
      size_t n = t.size();
      xt::xtensor<double, 2> data({n,3});
      for (int j = 0; j < n; ++j) {
        data(j, 0) = t[j].x;
        data(j, 1) = t[j].y;
        data(j, 2) = t[j].z;
      }
      std::string name = "coordinates_" + std::to_string(i);
      write_dataset(file_id, name.c_str(), data);
    }
    file_close(file_id);
  }

  // Clear particle tracks
  tracks.clear();
}

} // namespace openmc
