#include "openmc/ncrystal_interface.h"

#include "openmc/error.h"
#include "openmc/material.h"
#include "openmc/random_lcg.h"

namespace openmc {

//==============================================================================
// NCrystalMat implementation
//==============================================================================

NCrystalMat::NCrystalMat(const std::string& cfg)
  : cfg_( cfg ), proc_( cfg.c_str() )
{
}

double NCrystalMat::xs(const Particle& p) const
{
  // Calculate scattering XS per atom with NCrystal, only once per material
  double neutron_state[4] = { p.E(), p.u().x, p.u().y, p.u().z };
  return proc_.cross_section( neutron_state );
}

void NCrystalMat::scatter(Particle& p) const
{
  // Scatter with NCrystal, using the OpenMC RNG stream:
  uint64_t* seed = p.current_seed();
  std::function<double()> rng = [&seed]() { return prn(seed); };
  double neutron_state[4] = { p.E(), p.u().x, p.u().y, p.u().z };
  proc_.scatter(rng,neutron_state);
  // Modify attributes of particle
  p.E() = neutron_state[0];
  Direction u_old {p.u()};
  p.u() = Direction(neutron_state[1],neutron_state[2],neutron_state[3]);
  p.mu() = u_old.dot(p.u());
  p.event_mt() = ELASTIC;
}

//==============================================================================
// Functions
//==============================================================================

void ncrystal_update_micro(double xs, NuclideMicroXS& micro)
{
  if (micro.thermal > 0 || micro.thermal_elastic > 0) {
    fatal_error("S(a,b) treatment and NCrystal are not compatible.");
  }
  // remove free atom cross section
  // and replace it by scattering cross section per atom from NCrystal
  micro.total = micro.total - micro.elastic + xs;
  micro.elastic = xs;
}

} // namespace openmc
