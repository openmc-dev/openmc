#ifndef OPENMC_TALLIES_FILTER_ENERGY_H
#define OPENMC_TALLIES_FILTER_ENERGY_H

#include <vector>

#include <gsl/gsl>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Bins the incident neutron energy.
//==============================================================================

class EnergyFilter : public Filter
{
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~EnergyFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "energy";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const std::vector<double>& bins() const { return bins_; }
  void set_bins(gsl::span<const double> bins);

  bool matches_transport_groups() const { return matches_transport_groups_; }

protected:
  //----------------------------------------------------------------------------
  // Data members

  std::vector<double> bins_;

  //! True if transport group number can be used directly to get bin number
  bool matches_transport_groups_ {false};
};

//==============================================================================
//! Bins the outgoing neutron energy.
//!
//! Only scattering events use the get_all_bins functionality.  Nu-fission
//! tallies manually iterate over the filter bins.
//==============================================================================

class EnergyoutFilter : public EnergyFilter
{
public:
  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "energyout";}

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  std::string text_label(int bin) const override;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_ENERGY_H
