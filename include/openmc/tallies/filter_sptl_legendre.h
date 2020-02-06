#ifndef OPENMC_TALLIES_FILTER_SPTL_LEGENDRE_H
#define OPENMC_TALLIES_FILTER_SPTL_LEGENDRE_H

#include <string>

#include "openmc/tallies/filter.h"

namespace openmc {

enum class LegendreAxis {
  x, y, z
};

//==============================================================================
//! Gives Legendre moments of the particle's normalized position along an axis
//==============================================================================

class SpatialLegendreFilter : public Filter
{
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~SpatialLegendreFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "spatiallegendre";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  int order() const { return order_; }
  void set_order(int order);

  LegendreAxis axis() const { return axis_; }
  void set_axis(LegendreAxis axis);

  double min() const { return min_; }
  double max() const { return max_; }
  void set_minmax(double min, double max);

private:
  //----------------------------------------------------------------------------
  // Data members

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
