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
  // Methods
  void update_neutron_source(double k_eff) override;
  void normalize_scalar_flux_and_volumes(
    double total_active_distance_per_iteration) override;

  void batch_reset() override;
  void accumulate_iteration_flux() override;
  void output_to_vtk() const;
  double evaluate_flux_at_point(Position r, int64_t sr, int g) const override;

protected:
  //----------------------------------------------------------------------------
  // Methods
  void set_flux_to_flux_plus_source(int64_t sr, double volume, int g) override;
  void set_flux_to_zero(int64_t sr, int g) override;
  void set_flux_to_old_flux(int64_t sr, int g) override;

}; // class LinearSourceDomain

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_LINEAR_SOURCE_DOMAIN_H
