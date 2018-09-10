#include "openmc/distribution_multi.h"

#include <algorithm> // for move
#include <cmath>     // for sqrt, sin, cos, max

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/random_lcg.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// UnitSphereDistribution implementation
//==============================================================================

UnitSphereDistribution::UnitSphereDistribution(pugi::xml_node node)
{
  // Read reference directional unit vector
  if (check_for_node(node, "reference_uvw")) {
    auto u_ref = get_node_array<double>(node, "reference_uvw");
    if (u_ref.size() != 3)
      fatal_error("Angular distribution reference direction must have "
                  "three parameters specified.");
    u_ref_ = Direction(u_ref.data());
  }
}


//==============================================================================
// PolarAzimuthal implementation
//==============================================================================

PolarAzimuthal::PolarAzimuthal(Direction u, UPtrDist mu, UPtrDist phi) :
  UnitSphereDistribution{u}, mu_{std::move(mu)}, phi_{std::move(phi)} { }

PolarAzimuthal::PolarAzimuthal(pugi::xml_node node)
  : UnitSphereDistribution{node}
{
  if (check_for_node(node, "mu")) {
    pugi::xml_node node_dist = node.child("mu");
    mu_ = distribution_from_xml(node_dist);
  } else {
    mu_ = UPtrDist{new Uniform(-1., 1.)};
  }

  if (check_for_node(node, "phi")) {
    pugi::xml_node node_dist = node.child("phi");
    phi_ = distribution_from_xml(node_dist);
  } else {
    phi_ = UPtrDist{new Uniform(0.0, 2.0*PI)};
  }
}

Direction PolarAzimuthal::sample() const
{
  // Sample cosine of polar angle
  double mu = mu_->sample();
  if (mu == 1.0) return u_ref_;

  // Sample azimuthal angle
  double phi = phi_->sample();
  return rotate_angle(u_ref_, mu, &phi);
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
  return u_ref_;
}

} // namespace openmc
