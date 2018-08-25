#include "openmc/distribution_multi.h"

#include <algorithm> // for move
#include <cmath>     // for sqrt, sin, cos, max

#include "openmc/constants.h"
#include "openmc/math_functions.h"
#include "openmc/random_lcg.h"

namespace openmc {

//==============================================================================
// PolarAzimuthal implementation
//==============================================================================

PolarAzimuthal::PolarAzimuthal(Direction u, UPtrDist mu, UPtrDist phi) :
  UnitSphereDistribution{u}, mu_{std::move(mu)}, phi_{std::move(phi)} { }

Direction PolarAzimuthal::sample() const
{
  // Sample cosine of polar angle
  double mu = mu_->sample();
  if (mu == 1.0) return u_ref;

  // Sample azimuthal angle
  double phi = phi_->sample();
  return rotate_angle(u_ref, mu, &phi);
}

//==============================================================================
// Isotropic implementation
//==============================================================================

Direction Isotropic::sample() const
{
  double phi = 2.0*PI*prn();
  double mu = 2.0*prn() - 1.0;
  return {mu, std::sqrt(1.0 - mu*mu) * std::cos(phi),
      std::sqrt(1.0 - mu*mu) * std::sin(phi)};
}

//==============================================================================
// Monodirectional implementation
//==============================================================================

Direction Monodirectional::sample() const
{
  return u_ref;
}

} // namespace openmc
