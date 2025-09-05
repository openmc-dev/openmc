#ifndef OPENMC_TALLIES_FILTER_PARTICLE_H
#define OPENMC_TALLIES_FILTER_PARTICLE_H

#include "openmc/particle.h"
#include "openmc/span.h"
#include "openmc/tallies/filter.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Bins by type of particle (e.g. neutron, photon).
//==============================================================================

class ParticleFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~ParticleFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "particle"; }
  FilterType type() const override { return FilterType::PARTICLE; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const vector<ParticleType>& particles() const { return particles_; }

  void set_particles(span<ParticleType> particles);

private:
  //----------------------------------------------------------------------------
  // Data members

  vector<ParticleType> particles_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_PARTICLE_H
