#ifndef OPENMC_TALLIES_FILTER_MU_SURFACE_H
#define OPENMC_TALLIES_FILTER_MU_SURFACE_H

#include <gsl/gsl-lite.hpp>

#include "openmc/tallies/filter.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Bins the incoming-outgoing direction cosine.  This is only used for scatter
//! reactions.
//==============================================================================

class MuSurfaceFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~MuSurfaceFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "musurface"; }
  FilterType type() const override { return FilterType::MUSURFACE; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  void set_bins(gsl::span<double> bins);

private:
  //----------------------------------------------------------------------------
  // Data members

  vector<double> bins_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_MU_SURFACE_H