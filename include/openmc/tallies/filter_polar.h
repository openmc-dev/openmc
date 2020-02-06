#ifndef OPENMC_TALLIES_FILTER_POLAR_H
#define OPENMC_TALLIES_FILTER_POLAR_H

#include <cmath>
#include <vector>

#include <gsl/gsl>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Bins the incident neutron polar angle (relative to the global z-axis).
//==============================================================================

class PolarFilter : public Filter
{
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~PolarFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "polar";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  void set_bins(gsl::span<double> bins);

private:
  //----------------------------------------------------------------------------
  // Data members

  std::vector<double> bins_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_POLAR_H
