#ifndef OPENMC_TALLIES_FILTER_ENERGYFUNC_H
#define OPENMC_TALLIES_FILTER_ENERGYFUNC_H

#include "openmc/constants.h"
#include "openmc/span.h"
#include "openmc/tallies/filter.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Multiplies tally scores by an arbitrary function of incident energy
//! described by a piecewise linear-linear interpolation.
//==============================================================================

class EnergyFunctionFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  EnergyFunctionFilter() : Filter {} { n_bins_ = 1; }

  ~EnergyFunctionFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "energyfunction"; }
  FilterType type() const override { return FilterType::ENERGY_FUNCTION; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const vector<double>& energy() const { return energy_; }
  const vector<double>& y() const { return y_; }
  Interpolation interpolation() const { return interpolation_; }
  void set_data(span<const double> energy, span<const double> y);
  void set_interpolation(const std::string& interpolation);

private:
  //----------------------------------------------------------------------------
  // Data members

  //! Incident neutron energy interpolation grid.
  vector<double> energy_;

  //! Interpolant values.
  vector<double> y_;

  //! Interpolation scheme
  Interpolation interpolation_ {Interpolation::lin_lin};
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_ENERGYFUNC_H
