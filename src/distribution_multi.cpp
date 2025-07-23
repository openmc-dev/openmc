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

std::pair<Direction, double> PolarAzimuthal::sample(uint64_t* seed) const
{
  // Sample cosine of polar angle
  auto [mu, mu_wgt] = mu_->sample(seed);

  // Sample azimuthal angle
  auto [phi, phi_wgt] = phi_->sample(seed);
  if (mu == 1.0)
    return {u_ref_, mu_wgt * phi_wgt};

  return {rotate_angle(u_ref_, mu, &phi, seed), mu_wgt * phi_wgt};
}

std::pair<Direction, double> PolarAzimuthal::sample_as_bias(
  uint64_t* seed) const
{
  // Sample cosine of polar angle
  auto [mu, mu_wgt] = mu_->sample(seed);

  // Sample azimuthal angle
  auto [phi, phi_wgt] = phi_->sample(seed);

  double pdf_evaluation = mu_->evaluate(mu) * phi_->evaluate(phi);

  if (mu == 1.0)
    return {u_ref_, pdf_evaluation};

  return {rotate_angle(u_ref_, mu, &phi, seed), pdf_evaluation};
}

//==============================================================================
// Isotropic implementation
//==============================================================================

Isotropic::Isotropic(pugi::xml_node node) : UnitSphereDistribution {node}
{
  if (check_for_node(node, "bias")) {
    pugi::xml_node bias_node = node.child("bias");
    std::string bias_type = get_node_value(bias_node, "type", true, true);
    if (bias_type != "mu-phi") {
      openmc::fatal_error(
        "Isotropic distributions may only be biased by a PolarAzimuthal.");
    }
    auto bias = std::make_unique<PolarAzimuthal>(bias_node);
    if (bias->mu()->bias() || bias->phi()->bias()) {
      openmc::fatal_error(
        "Attempted to bias Isotropic distribution with a biased PolarAzimuthal "
        "distribution. Please ensure bias distributions are unbiased.");
    }
    this->set_bias(std::move(bias));
  }
}

Direction isotropic_direction(uint64_t* seed)
{
  double phi = uniform_distribution(0., 2.0 * PI, seed);
  double mu = uniform_distribution(-1., 1., seed);
  return {mu, std::sqrt(1.0 - mu * mu) * std::cos(phi),
    std::sqrt(1.0 - mu * mu) * std::sin(phi)};
}

std::pair<Direction, double> Isotropic::sample(uint64_t* seed) const
{
  if (bias()) {
    auto [val, eval] = bias()->sample_as_bias(seed);
    return {val, 1.0 / (4.0 * PI * eval)};
  } else {
    return {isotropic_direction(seed), 1.0};
  }
}

//==============================================================================
// Monodirectional implementation
//==============================================================================

std::pair<Direction, double> Monodirectional::sample(uint64_t* seed) const
{
  return {u_ref_, 1.0};
}

} // namespace openmc
