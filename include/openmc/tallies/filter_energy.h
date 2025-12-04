#ifndef OPENMC_TALLIES_FILTER_ENERGY_H
#define OPENMC_TALLIES_FILTER_ENERGY_H

#include "openmc/particle.h"
#include "openmc/span.h"
#include "openmc/tallies/filter.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Bins the incident neutron energy.
//==============================================================================

class EnergyFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~EnergyFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "energy"; }
  FilterType type() const override { return FilterType::ENERGY; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const vector<double>& bins() const { return bins_; }
  void set_bins(span<const double> bins);

  bool matches_transport_groups() const { return matches_transport_groups_; }

protected:
  //----------------------------------------------------------------------------
  // Data members

  vector<double> bins_;

  //! True if transport group number can be used directly to get bin number
  bool matches_transport_groups_ {false};
};

//==============================================================================
//! Bins the outgoing neutron energy.
//!
//! Only scattering events use the get_all_bins functionality.  Nu-fission
//! tallies manually iterate over the filter bins.
//==============================================================================

class EnergyoutFilter : public EnergyFilter {
public:
  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "energyout"; }
  FilterType type() const override { return FilterType::ENERGY_OUT; }

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  std::string text_label(int bin) const override;
};

//==============================================================================
//! Bins the outgoing energy of secondary particles
//!
//! This is used to get the photon production matrix for multigroup photon
//! transport. It could also be used to find the energy distribution of
//! neutron secondaries or others, for example.
//!
//! Using anything other than analog estimators here would be complicated
//==============================================================================

class SecondaryEnergyFilter : public EnergyFilter {
public:
  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "secondaryenergy"; }
  FilterType type() const override { return FilterType::SECONDARY_ENERGY; }

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  std::string text_label(int bin) const override;

  void from_xml(pugi::xml_node node) override;

protected:
  // This filter could simultaneously use different particle types, for
  // example creating one set of energy bins for photons and another for
  // electrons. However, for typical use only photons are of interest.
  // Covering multiple particles can be done by defining separate tallies.
  ParticleType secondary_type_ {ParticleType::photon};
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_ENERGY_H
