#ifndef OPENMC_TALLIES_FILTER_PARTICLE_H
#define OPENMC_TALLIES_FILTER_PARTICLE_H

#include <vector>

#include "openmc/particle.h"
#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Bins by type of particle (e.g. neutron, photon).
//==============================================================================

class ParticleFilter : public Filter
{
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~ParticleFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "particle";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const std::vector<Particle::Type>& particles() const { return particles_; }

  void set_particles(gsl::span<Particle::Type> particles);

private:
  //----------------------------------------------------------------------------
  // Data members

  std::vector<Particle::Type> particles_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_PARTICLE_H
