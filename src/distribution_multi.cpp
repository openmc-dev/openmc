#include "openmc/distribution_multi.h"

#include <algorithm> // for move
#include <cmath>     // for sqrt, sin, cos, max

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/random_dist.h"
#include "openmc/random_lcg.h"
#include "openmc/xml_interface.h"

namespace openmc {

unique_ptr<UnitSphereDistribution> UnitSphereDistribution::create(
  pugi::xml_node node)
{
  // Check for type of angular distribution
  std::string type;
  if (check_for_node(node, "type"))
    type = get_node_value(node, "type", true, true);
  if (type == "isotropic") {
    return UPtrAngle {new Isotropic()};
  } else if (type == "monodirectional") {
    return UPtrAngle {new Monodirectional(node)};
  } else if (type == "mu-phi") {
    return UPtrAngle {new PolarAzimuthal(node)};
  } else {
    fatal_error(fmt::format(
      "Invalid angular distribution for external source: {}", type));
  }
}

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

PolarAzimuthal::PolarAzimuthal(Direction u, UPtrDist mu, UPtrDist phi)
  : UnitSphereDistribution {u}, mu_ {std::move(mu)}, phi_ {std::move(phi)}
{}

PolarAzimuthal::PolarAzimuthal(pugi::xml_node node)
  : UnitSphereDistribution {node}
{
  // Read reference directional unit vector
  if (check_for_node(node, "reference_vwu")) {
    auto v_ref = get_node_array<double>(node, "reference_vwu");
    if (v_ref.size() != 3)
      fatal_error("Angular distribution reference v direction must have "
                  "three parameters specified.");
    v_ref_ = Direction(v_ref.data());
  }
  w_ref_ = u_ref_.cross(v_ref_);
  if (check_for_node(node, "mu")) {
    pugi::xml_node node_dist = node.child("mu");
    mu_ = distribution_from_xml(node_dist);
  } else {
    mu_ = UPtrDist {new Uniform(-1., 1.)};
  }

  if (check_for_node(node, "phi")) {
    pugi::xml_node node_dist = node.child("phi");
    phi_ = distribution_from_xml(node_dist);
  } else {
    phi_ = UPtrDist {new Uniform(0.0, 2.0 * PI)};
  }
}

Direction PolarAzimuthal::sample(uint64_t* seed) const
{
  // Sample cosine of polar angle
  double mu = mu_->sample(seed);
  if (mu == 1.0)
    return u_ref_;
  if (mu == -1.0)
    return -u_ref_;
    
  // Sample azimuthal angle
  double phi = phi_->sample(seed);
  
  double f = std::sqrt(1-mu*mu);
  
  return mu*u_ref_ + f*std::cos(phi)*v_ref_ + f*std::sin(phi)*w_ref_;
}

//==============================================================================
// Isotropic implementation
//==============================================================================

Direction isotropic_direction(uint64_t* seed)
{
  double phi = uniform_distribution(0., 2.0 * PI, seed);
  double mu = uniform_distribution(-1., 1., seed);
  return {mu, std::sqrt(1.0 - mu * mu) * std::cos(phi),
    std::sqrt(1.0 - mu * mu) * std::sin(phi)};
}

Direction Isotropic::sample(uint64_t* seed) const
{
  return isotropic_direction(seed);
}

//==============================================================================
// Monodirectional implementation
//==============================================================================

Direction Monodirectional::sample(uint64_t* seed) const
{
  return u_ref_;
}

} // namespace openmc
