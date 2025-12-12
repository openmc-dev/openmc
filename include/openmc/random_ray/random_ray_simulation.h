#ifndef OPENMC_RANDOM_RAY_SIMULATION_H
#define OPENMC_RANDOM_RAY_SIMULATION_H

#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/linear_source_domain.h"

namespace openmc {

/*
 * The RandomRaySimulation class encompasses data and methods for running a
 * random ray simulation.
 */

class RandomRaySimulation {
public:
  //----------------------------------------------------------------------------
  // Constructors
  RandomRaySimulation();

  //----------------------------------------------------------------------------
  // Methods
  void compute_segment_correction_factors();
  void apply_fixed_sources_and_mesh_domains();
  void prepare_fixed_sources_adjoint();
  void simulate();
  void output_simulation_results() const;
  void instability_check(
    int64_t n_hits, double k_eff, double& avg_miss_rate) const;
  void print_results_random_ray() const;
  int64_t total_geometric_intersections() const
  {
    return total_geometric_intersections_;
  }

  //----------------------------------------------------------------------------
  // Accessors
  FlatSourceDomain* domain() const { return domain_.get(); }

private:
  //----------------------------------------------------------------------------
  // Data members

  // Contains all flat source region data
  unique_ptr<FlatSourceDomain> domain_;

  // Tracks the average FSR miss rate for analysis and reporting
  double avg_miss_rate_ {0.0};

  // Tracks the total number of geometric intersections by all rays for
  // reporting
  uint64_t total_geometric_intersections_ {0};

  // Number of energy groups
  int negroups_;

}; // class RandomRaySimulation

//============================================================================
//! Non-member functions
//============================================================================

void openmc_run_random_ray();
void validate_random_ray_inputs();
void openmc_reset_random_ray();

//! Write data related to randaom ray to statepoint
//! \param[in] group HDF5 group
void write_random_ray_hdf5(hid_t group);

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_SIMULATION_H
