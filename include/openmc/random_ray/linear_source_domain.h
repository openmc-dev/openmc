#ifndef OPENMC_RANDOM_RAY_LINEAR_SOURCE_DOMAIN_H
#define OPENMC_RANDOM_RAY_LINEAR_SOURCE_DOMAIN_H

#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/moment_matrix.h"

#include "openmc/openmp_interface.h"
#include "openmc/position.h"
#include "openmc/source.h"

namespace openmc {

/*
 * The LinearSourceDomain class encompasses data and methods for storing
 * scalar flux and source region for all linear source regions in a
 * random ray simulation domain.
 */

class LinearSourceDomain : public FlatSourceDomain {
public:
  //----------------------------------------------------------------------------
  // Constructors
  LinearSourceDomain();

  //----------------------------------------------------------------------------
  // Methods
  void update_neutron_source(double k_eff) override;
  double compute_k_eff(double k_eff_old) const;
  void normalize_scalar_flux_and_volumes(
    double total_active_distance_per_iteration) override;
  int64_t add_source_to_scalar_flux() override;
  void batch_reset() override;
  void convert_source_regions_to_tallies();
  void reset_tally_volumes();
  void random_ray_tally();
  void accumulate_iteration_flux() override;
  void output_to_vtk() const;
  void all_reduce_replicated_source_regions() override;
  void convert_external_sources();
  void count_external_source_regions();
  void flux_swap() override;
  double evaluate_flux_at_point(Position r, int64_t sr, int g) const override;

  void compute_first_collided_external_source() override;
  void compute_uncollided_scalar_flux() override;
  void compute_first_collided_flux() override;
  void normalize_uncollided_scalar_flux(double number_of_particles) override;
  void update_volume_uncollided_flux() override;
  void uncollided_moments() override;
  void batch_reset_fc() override;

  //----------------------------------------------------------------------------
  // Public Data members

  vector<MomentArray> source_gradients_;
  vector<MomentArray> flux_moments_old_;
  vector<MomentArray> flux_moments_new_;
  vector<MomentArray> flux_moments_t_;
  // First Collided Method
  vector<MomentArray> flux_moments_uncollided_;
  vector<MomentArray> flux_moments_first_collided_;
  vector<MomentArray> external_source_gradients_;


  vector<Position> centroid_;
  vector<Position> centroid_iteration_;
  vector<Position> centroid_t_;
  vector<MomentMatrix> mom_matrix_;
  vector<MomentMatrix> mom_matrix_t_;

}; // class LinearSourceDomain

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_LINEAR_SOURCE_DOMAIN_H
