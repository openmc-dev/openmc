#ifndef OPENMC_RANDOM_RAY_LINEAR_SOURCE_DOMAIN_H
#define OPENMC_RANDOM_RAY_LINEAR_SOURCE_DOMAIN_H

#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/symmetric_matrix.h"

#include "openmc/openmp_interface.h"
#include "openmc/position.h"
#include "openmc/source.h"

namespace openmc {

/*
 * The LinearSourceDomain class encompasses data and methods for storing
 * scalar flux and source region for all flat source regions in a
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
  double compute_k_eff(double k_eff_old) const; // override
  void normalize_scalar_flux_and_volumes(
    double total_active_distance_per_iteration) override;
  int64_t add_source_to_scalar_flux() override;
  void batch_reset() override;
  void convert_source_regions_to_tallies();
  void reset_tally_volumes();
  void random_ray_tally();
  void accumulate_iteration_flux() override;
  void output_to_vtk() const;
  void all_reduce_replicated_source_regions(); //override
  void convert_external_sources();
  void count_external_source_regions();
  void flux_swap() override;
  double evaluate_flux_at_point(const Position r, const int64_t sr, const int g) const override;


  //----------------------------------------------------------------------------
  // Public Data members

  vector<float> flux_x_new_;
  vector<float> flux_x_old_;
    vector<float> flux_x_t_;

  vector<float> source_x_;
  vector<float> flux_y_new_;
    vector<float> flux_y_t_;
  vector<float> flux_y_old_;
  vector<float> source_y_;
  vector<float> flux_z_new_;
  vector<float> flux_z_t_;
  vector<float> flux_z_old_;
  vector<float> source_z_;
  vector<double> centroid_;
  vector<double> centroid_t_;
  vector<SymmetricMatrix> mom_matrix_;
  vector<SymmetricMatrix> mom_matrix_t_;

private:
  //----------------------------------------------------------------------------
  // Methods

  // ...

  //----------------------------------------------------------------------------
  // Private data members
  
  // ...

}; // class LinearSourceDomain

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_LINEAR_SOURCE_DOMAIN_H