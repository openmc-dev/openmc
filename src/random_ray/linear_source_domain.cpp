#include "openmc/random_ray/linear_source_domain.h"

#include "openmc/cell.h"
#include "openmc/geometry.h"
#include "openmc/material.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/output.h"
#include "openmc/plot.h"
#include "openmc/random_ray/random_ray.h"
#include "openmc/simulation.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/timer.h"

#include <cstdio>

namespace openmc {

//==============================================================================
// LinearSourceDomain implementation
//==============================================================================

LinearSourceDomain::LinearSourceDomain() : FlatSourceDomain()
{
    // The base class constructor is already called, so just need to initialize
    // the extra linear source members below
    // ...
}

void LinearSourceDomain::batch_reset()
{
    // We will first call the base class method as it needs to reset these
    // variables as well
    FlatSourceDomain::batch_reset();
    // And then we add in the stuff specific to the linear source class as well
    // ...
}

void LinearSourceDomain::accumulate_iteration_flux()
{
    // Fully reimplement this function for linear source
    // ...
}

void LinearSourceDomain::update_neutron_source(double k_eff)
{
    // Fully reimplement this function for linear source
    // ...
}

void LinearSourceDomain::normalize_scalar_flux_and_volumes(
  double total_active_distance_per_iteration)
  {
    // Fully reimplement this function for linear source
    // ...
  }

  int64_t LinearSourceDomain::add_source_to_scalar_flux()
{
    int64_t n_hits = 0;
    // Fully reimplement this function for linear source
    // ...
    return n_hits;
}

double LinearSourceDomain::compute_k_eff(double k_eff_old) const
{
    // Fully reimplement this function for linear source
    // ...
    return 1.0;
}

//void LinearSourceDomain::random_ray_tally()
//{
    // Can we just use the base class method? (Do we need this method at all?)
    // ...
//}

void LinearSourceDomain::all_reduce_replicated_source_regions()
{
    // We will first call the base class method as it needs to reduce these
    // variables as well
    FlatSourceDomain::all_reduce_replicated_source_regions();
    // Then we add in the stuff specific to the linear source class
    // ...
}

//oid LinearSourceDomain::output_to_vtk() const
//{
    // Can we just use the base class method? (Do we need this method at all?)
    // ...
//}

//void LinearSourceDomain::apply_external_source_to_source_region(
//  Discrete* discrete, double strength_factor, int64_t source_region)
//{
    // Can we just use the base class method? (Do we need this method at all?)
    // ...
//}

} // namespace openmc