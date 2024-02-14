#ifndef OPENMC_RANDOM_RAY_ITERATION_H
#define OPENMC_RANDOM_RAY_ITERATION_H

#include "openmc/random_ray/random_ray.h"
#include "openmc/random_ray/source_region.h"

namespace openmc {

/*
 * The RandomRaySimulation class encompasses data and methods for running a
 * random ray simulation.
 */
  
class RandomRaySimulation {
public:
  //==========================================================================
  // Constructors
  RandomRaySimulation();
  
  //==========================================================================
  // Methods
  
  void simulate();

  void all_reduce_random_ray_batch_results(bool mapped_all_tallies);
  void find_random_ray_sampling_source_index();
  void reduce_simulation_statistics();
  void output_simulation_results();
  void instability_check(int64_t n_hits, double k_eff, double& avg_miss_rate);

  //==========================================================================
  // Data

  FlatSourceDomain domain_;

  // Random ray eigenvalue
  double k_eff_ {1.0};
  
  // Tracks the average FSR miss rate for analysis and reporting
  double avg_miss_rate_ {0.0};

  // Tracks the total number of geometric intersections by all rays for
  // reporting
  uint64_t total_geometric_intersections_ {0};

  int sampling_source_;

  bool mapped_all_tallies_ {false};

  int negroups_;

}; // class RandomRaySimulation


//==========================================================================
// Non-Method Functions

void openmc_run_random_ray();
void validate_random_ray_inputs();


} // namespace openmc

#endif // OPENMC_RANDOM_RAY_ITERATION_H
