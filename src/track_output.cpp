#include "openmc/track_output.h"

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/position.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"

#include <fmt/core.h>
#include "xtensor/xtensor.hpp"

#include <cstddef> // for size_t
#include <string>
#include <vector>

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

//==============================================================================
// Non-member functions
//==============================================================================

void add_particle_track(Particle& p)
{
  p.tracks_.emplace_back();
}

void write_particle_track(Particle& p)
{
  p.tracks_.back().push_back(p.r());
}

void finalize_particle_track(Particle& p)
{
  std::string filename = fmt::format("{}track_{}_{}_{}.h5",
    settings::path_output, simulation::current_batch, simulation::current_gen,
    p.id_);

  // Determine number of coordinates for each particle
  std::vector<int> n_coords;
  for (auto& coords : p.tracks_) {
    n_coords.push_back(coords.size());
  }

  #pragma omp critical (FinalizeParticleTrack)
  {
    hid_t file_id = file_open(filename, 'w');
    write_attribute(file_id, "filetype", "track");
    write_attribute(file_id, "version", VERSION_TRACK);
    write_attribute(file_id, "n_particles", p.tracks_.size());
    write_attribute(file_id, "n_coords", n_coords);
    for (auto i = 1; i <= p.tracks_.size(); ++i) {
      const auto& t {p.tracks_[i-1]};
      size_t n = t.size();
      xt::xtensor<double, 2> data({n,3});
      for (int j = 0; j < n; ++j) {
        data(j, 0) = t[j].x;
        data(j, 1) = t[j].y;
        data(j, 2) = t[j].z;
      }
      std::string name = fmt::format("coordinates_{}", i);
      write_dataset(file_id, name.c_str(), data);
    }
    file_close(file_id);
  }

  // Clear particle tracks
  p.tracks_.clear();
}

} // namespace openmc
