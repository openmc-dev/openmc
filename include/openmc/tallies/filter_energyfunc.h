#ifndef OPENMC_TALLIES_FILTER_ENERGYFUNC_H
#define OPENMC_TALLIES_FILTER_ENERGYFUNC_H

#include <vector>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Multiplies tally scores by an arbitrary function of incident energy
//! described by a piecewise linear-linear interpolation.
//==============================================================================

class EnergyFunctionFilter : public Filter
{
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  EnergyFunctionFilter()
    : Filter {}
  {
    n_bins_ = 1;
  }

  ~EnergyFunctionFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "energyfunction";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const std::vector<double>& energy() const { return energy_; }
  const std::vector<double>& y() const { return y_; }
  void set_data(gsl::span<const double> energy, gsl::span<const double> y);

private:
  //----------------------------------------------------------------------------
  // Data members

  //! Incident neutron energy interpolation grid.
  std::vector<double> energy_;

  //! Interpolant values.
  std::vector<double> y_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_ENERGYFUNC_H
