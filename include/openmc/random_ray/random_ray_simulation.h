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
  void apply_fixed_sources_and_mesh_domains();
  void prepare_fixed_sources_adjoint();
  void print_random_ray_headers();
  void run_single_simulation();
  void random_ray_adjoint();
  void kinetic_initial_condition();
  void kinetic_single_time_step(int i);
  void simulate();
  void output_simulation_results() const;
  void instability_check(
    int64_t n_hits, double k_eff, double& avg_miss_rate) const;
  void print_results_random_ray() const;

  //----------------------------------------------------------------------------
  // Accessors
  FlatSourceDomain* domain() const { return domain_.get(); }

private:
  //----------------------------------------------------------------------------
  // Data Members

  // Contains all flat source region data
  unique_ptr<FlatSourceDomain> domain_;

  // Tracks the average FSR miss rate for analysis and reporting
  double avg_miss_rate_ {0.0};

  // Tracks the total number of geometric intersections by all rays for
  // reporting
  uint64_t total_geometric_intersections_ {0};

  // Number of energy groups
  int negroups_;

  // Number of delay groups
  int ndgroups_;

  // Toggle for first simulation
  bool is_first_simulation_;

  // Flag for adjoint simulation;
  bool adjoint_needed_;

  //----------------------------------------------------------------------------
  // Data Members for kinetic simulations

  double static_avg_k_eff_;
  vector<double> static_k_eff_;
  vector<double> static_fission_rate_;

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
void print_random_ray_headers(bool& adjoint_needed);

// Functions for kinetic simulations
void set_time_dependent_settings();
void rename_time_step_file(
  std::string base_filename, std::string extension, int i);

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_SIMULATION_H
