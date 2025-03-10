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
  void prepare_fixed_sources_adjoint(vector<double>& forward_flux,
    SourceRegionContainer& forward_source_regions,
    SourceRegionContainer& forward_base_source_regions,
    std::unordered_map<SourceRegionKey, int64_t, SourceRegionKey::HashFunctor>&
      forward_source_region_map);
  void simulate();
  void output_simulation_results() const;
  void instability_check(
    int64_t n_hits, double k_eff, double& avg_miss_rate) const;
  void print_results_random_ray(uint64_t total_geometric_intersections,
    double avg_miss_rate, int negroups, int64_t n_source_regions,
    int64_t n_external_source_regions) const;

  //----------------------------------------------------------------------------
  // Accessors
  FlatSourceDomain* domain() const { return domain_.get(); }

private:
  //----------------------------------------------------------------------------
  // Data members

  // Contains all flat source region data
  unique_ptr<FlatSourceDomain> domain_;

  // Random ray eigenvalue
  double k_eff_ {1.0};

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

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_SIMULATION_H
