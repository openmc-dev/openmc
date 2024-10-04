#ifndef OPENMC_TALLIES_FILTER_MU_SURFACE_H
#define OPENMC_TALLIES_FILTER_MU_SURFACE_H

#include <gsl/gsl-lite.hpp>

#include "openmc/tallies/filter_mu.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Bins the incoming-outgoing direction cosine.  This is only used for surface
//! crossings.
//==============================================================================

class MuSurfaceFilter : public MuFilter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~MuSurfaceFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "musurface"; }
  FilterType type() const override { return FilterType::MUSURFACE; }

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_MU_SURFACE_H
