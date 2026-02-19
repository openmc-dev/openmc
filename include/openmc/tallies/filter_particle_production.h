#ifndef OPENMC_TALLIES_FILTER_PARTICLE_PRODUCTION_H
#define OPENMC_TALLIES_FILTER_PARTICLE_PRODUCTION_H

#include <unordered_map>

#include "openmc/particle.h"
#include "openmc/tallies/filter.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Bins the outgoing energy of secondary particles
//!
//! This filter bins secondary particles by type and, optionally, by energy. It
//! can be used to get the photon production matrix for multigroup photon
//! transport, the energy distribution of secondary neutrons, etc. The weight
//! that is applied is equal to the weight of the secondary particle. Thus, to
//! get secondary production it should be used in conjunction with the "events"
//! score.
//==============================================================================

class ParticleProductionFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "particleproduction"; }
  FilterType type() const override { return FilterType::PARTICLE_PRODUCTION; }

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  std::string text_label(int bin) const override;

  void from_xml(pugi::xml_node node) override;

  void to_statepoint(hid_t filter_group) const override;

private:
  //----------------------------------------------------------------------------
  // Data members

  vector<ParticleType> secondary_types_; //!< Types of secondary particles

  //! Map from PDG number to index in secondary_types_ for O(1) lookup
  std::unordered_map<int32_t, int> type_to_index_;

  vector<double> energy_bins_; //!< Energy bin boundaries (optional)
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_PARTICLE_PRODUCTION_H
