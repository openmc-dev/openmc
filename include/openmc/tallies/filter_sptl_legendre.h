#ifndef OPENMC_TALLIES_FILTER_SPTL_LEGENDRE_H
#define OPENMC_TALLIES_FILTER_SPTL_LEGENDRE_H

#include <string>

#include "openmc/tallies/filter.h"

namespace openmc {

//TODO: those integer values are not needed when Fortran interop is removed
enum class LegendreAxis {
  x = 1, y = 2, z = 3
};

//==============================================================================
//! Gives Legendre moments of the particle's normalized position along an axis
//==============================================================================

class SpatialLegendreFilter : public Filter
{
public:
  ~SpatialLegendreFilter() = default;

  std::string type() const override {return "spatiallegendre";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  int order_;

  //! The Cartesian coordinate axis that the Legendre expansion is applied to.
  LegendreAxis axis_;

  //! The minimum coordinate along the reference axis that the expansion covers.
  double min_;

  //! The maximum coordinate along the reference axis that the expansion covers.
  double max_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_SPTL_LEGENDRE_H
