#ifndef OPENMC_RANDOM_RAY_SOURCE_REGION_H
#define OPENMC_RANDOM_RAY_SOURCE_REGION_H

#include <vector>

#include "openmc/openmp_interface.h"
#include "openmc/position.h"

namespace openmc {

/*
 * The FlatSourceDomain class encompasses data and methods for storing
 * scalar flux and source region for all flat source regions in a 
 * random ray simulation domain.
 */

class FlatSourceDomain {
public:
  //==========================================================================
  // Constructors
  FlatSourceDomain();
  
  
  //==========================================================================
  // Methods
  void update_neutron_source(double k_eff);
  double compute_k_eff(double k_eff_old);
  void normalize_scalar_flux_and_volumes();
  int64_t add_source_to_scalar_flux();
  double calculate_miss_rate();
  void batch_reset();
  void convert_source_regions_to_tallies();
  void random_ray_tally();
  void accumulate_iteration_flux();

  //==========================================================================
  // Data

  // Scalars
  int negroups_; // Number of energy groups in simulation
  int64_t n_source_elements_; // Total number of source regions in the model
                                // times the number of energy groups
  int64_t n_source_regions_;  // Total number of source regions in the model

  bool mapped_all_tallies_;

  std::vector<std::vector<TallyTask>> tally_task_;

  // 1D arrays representing values for each OpenMC "Cell"
  std::vector<int64_t> source_region_offsets_;

  // 1D arrays reprenting values for all source regions
  std::vector<OpenMPMutex> lock_;
  std::vector<int> material_;
  std::vector<int> position_recorded_;
  std::vector<Position> position_;
  std::vector<double> volume_;
  std::vector<double> volume_t_;
  std::vector<int> was_hit_;

  // 2D arrays stored in 1D representing values for all source regions x energy
  // groups
  std::vector<float> scalar_flux_new_;
  std::vector<float> scalar_flux_old_;
  std::vector<float> scalar_flux_final_;
  std::vector<float> source_;

}; // class FlatSourceDomain

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_SOURCE_REGION_H
