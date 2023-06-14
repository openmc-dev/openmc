#include "openmc/cell.h"
#include "openmc/random_ray/source_region.h"
#include "openmc/random_ray/tally_convert.h"
#include "openmc/mgxs_interface.h"
#include <cassert>

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace random_ray {

// Scalars
int64_t n_source_elements {0}; // Total number of source regions in the model times the number of energy groups
int64_t n_source_regions {0}; // Total number of source regions in the model

// 1D arrays representing values for each OpenMC "Cell"
std::vector<int64_t> source_region_offsets;

// 1D arrays reprenting values for all source regions
std::vector<OpenMPMutex> lock;
std::vector<int> material;
std::vector<int> position_recorded;
std::vector<Position> position;
std::vector<double> volume;
std::vector<double> volume_t;
std::vector<int> was_hit;

// 1D arrays representing values for all source regions x energy groups
std::vector<float> scalar_flux_new;
std::vector<float> scalar_flux_old;
std::vector<float> source;


} // namespace random_ray

//==============================================================================
// Non-method functions
//==============================================================================

void initialize_source_regions()
{
  int negroups = data::mg.num_energy_groups_;

  // Count the number of source regions, compute
  // the cell offset indices, and store the material type The reason for the offsets is that
  // some cell types may not have material fills, and therefore
  // do not produce FSRs. Thus, we cannot index into the global
  // arrays directly
  for (auto&& c : model::cells) {
    if (c->type_ != Fill::MATERIAL) {
      random_ray::source_region_offsets.push_back(-1);
    } else {
      random_ray::source_region_offsets.push_back(random_ray::n_source_regions);
      random_ray::n_source_regions += c->n_instances_;
      random_ray::n_source_elements += c->n_instances_ * negroups;
    }
  }

  // Initialize cell-wise arrays
  random_ray::lock.resize(random_ray::n_source_regions);
  random_ray::material.resize(random_ray::n_source_regions);
  random_ray::position_recorded.assign(random_ray::n_source_regions, 0);
  random_ray::position.resize(random_ray::n_source_regions);
  random_ray::volume.assign(random_ray::n_source_regions, 0.0);
  random_ray::volume_t.assign(random_ray::n_source_regions, 0.0);
  random_ray::was_hit.assign(random_ray::n_source_regions, 0);

  // Initialize element-wise arrays
  random_ray::scalar_flux_new.assign(random_ray::n_source_elements, 0.0);
  random_ray::scalar_flux_old.assign(random_ray::n_source_elements, 1.0);
  random_ray::source.resize(random_ray::n_source_elements);
  random_ray::tally_task.resize(random_ray::n_source_elements);

  // Initialize material array
  int64_t source_region_id = 0;
  for (int i =  0; i < model::cells.size(); i++) {
    Cell& cell = *model::cells[i];
    if (cell.type_ == Fill::MATERIAL) {
      int material = cell.material_[0];
      for (int j = 0; j < cell.n_instances_; j++) {
        random_ray::material[source_region_id++] = material;;
      }
    }
  }
  assert(source_region_id == random_ray::n_source_regions);
}

} // namespace openmc
