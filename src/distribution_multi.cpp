#include "openmc/distribution_multi.h"

#include <algorithm> // for move, clamp
#include <cmath>     // for sqrt, sin, cos, max

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/mesh.h"
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
    return UPtrAngle {new Isotropic(node)};
  } else if (type == "monodirectional") {
    return UPtrAngle {new Monodirectional(node)};
  } else if (type == "mu-phi") {
    return UPtrAngle {new PolarAzimuthal(node)};
  } else if (type == "mesh-angular") {
    return UPtrAngle {new MeshAngular(node)};
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
    u_ref_ /= u_ref_.norm();
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
    v_ref_ /= v_ref_.norm();
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

std::pair<Direction, double> PolarAzimuthal::sample(uint64_t* seed) const
{
  return sample_impl(seed, false);
}

std::pair<Direction, double> PolarAzimuthal::sample_as_bias(
  uint64_t* seed) const
{
  return sample_impl(seed, true);
}

std::pair<Direction, double> PolarAzimuthal::sample_impl(
  uint64_t* seed, bool return_pdf) const
{
  // Sample cosine of polar angle
  auto [mu, mu_wgt] = mu_->sample(seed);

  // Sample azimuthal angle
  auto [phi, phi_wgt] = phi_->sample(seed);

  // Compute either the PDF value or the importance weight
  double weight =
    return_pdf ? (mu_->evaluate(mu) * phi_->evaluate(phi)) : (mu_wgt * phi_wgt);

  if (mu == 1.0)
    return {u_ref_, weight};
  if (mu == -1.0)
    return {-u_ref_, weight};

  double f = std::sqrt(1 - mu * mu);
  return {mu * u_ref_ + f * std::cos(phi) * v_ref_ + f * std::sin(phi) * w_ref_,
    weight};
}

double PolarAzimuthal::evaluate(Direction u) const
{
  double mu = std::clamp(u.dot(u_ref_), -1.0, 1.0);
  double phi = 0.0;
  double sin_theta_sq = std::max(0.0, 1.0 - mu * mu);
  if (sin_theta_sq > 0.0) {
    double sin_theta = std::sqrt(sin_theta_sq);
    double cos_phi = u.dot(v_ref_) / sin_theta;
    double sin_phi = u.dot(w_ref_) / sin_theta;
    phi = std::atan2(sin_phi, cos_phi);
    if (phi < 0.0)
      phi += 2.0 * PI;
  }
  return mu_->evaluate(mu) * phi_->evaluate(phi);
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

double Isotropic::evaluate(Direction u) const
{
  return 1.0 / (4.0 * PI);
}

//==============================================================================
// Monodirectional implementation
//==============================================================================

std::pair<Direction, double> Monodirectional::sample(uint64_t* seed) const
{
  return {u_ref_, 1.0};
}

//==============================================================================
// MeshAngular implementation
//==============================================================================

MeshAngular::MeshAngular(pugi::xml_node node) : UnitSphereDistribution {node}
{
  auto spatial_type = get_node_value(node, "type", true, true);
  if (spatial_type != "mesh_angular") {
    fatal_error(
      fmt::format("Incorrect angular type '{}' for a MeshAngular distribution",
        spatial_type));
  }

  int32_t mesh_id = std::stoi(get_node_value(node, "mesh_id"));
  // Get pointer to spatial distribution
  mesh_idx_ = model::mesh_map.at(mesh_id);

  check_element_types();

  size_t n_bins = this->n_sources();
  std::vector<double> strengths(n_bins, 1.0);

  // Create cdfs for sampling for an element over a mesh
  if (check_for_node(node, "strengths")) {
    strengths = get_node_array<double>(node, "strengths");
    if (strengths.size() != n_bins) {
      fatal_error(
        fmt::format("Number of entries in the source strengths array {} does "
                    "not match the number of entities in mesh {} ({}).",
          strengths.size(), mesh_id, n_bins));
    }
  }

  elem_idx_dist_.assign(strengths);

  if (check_for_node(node, "bias")) {
    pugi::xml_node bias_node = node.child("bias");

    if (check_for_node(bias_node, "strengths")) {
      std::vector<double> bias_strengths(n_bins, 1.0);
      bias_strengths = get_node_array<double>(node, "strengths");

      if (bias_strengths.size() != n_bins) {
        fatal_error(
          fmt::format("Number of entries in the bias strengths array {} does "
                      "not match the number of entities in mesh {} ({}).",
            bias_strengths.size(), mesh_id, n_bins));
      }

      // Compute importance weights
      weight_ = compute_importance_weights(strengths, bias_strengths);

      // Re-initialize DiscreteIndex with bias strengths for sampling
      elem_idx_dist_.assign(bias_strengths);
    } else {
      fatal_error(fmt::format(
        "Bias node for mesh {} found without strengths array.", mesh_id));
    }
  }
}

void MeshAngular::check_element_types() const
{
  const auto trimesh_ptr =
    dynamic_cast<const UnitSphereTriangularMesh*>(this->mesh());
  if (!trimesh_ptr) {
    fatal_error("MeshAngular distribution must use a mesh of type "
                "UnitSphereTriangularMesh.");
  }
}

int32_t MeshAngular::sample_element_index(uint64_t* seed) const
{
  return elem_idx_dist_.sample(seed);
}

std::pair<int32_t, Direction> MeshAngular::sample_mesh(uint64_t* seed) const
{
  // Sample the CDF defined in initialization above
  int32_t elem_idx = this->sample_element_index(seed);
  return {elem_idx, mesh()->sample_element(elem_idx, seed)};
}

std::pair<Direction, double> MeshAngular::sample(uint64_t* seed) const
{
  auto [elem_idx, u] = this->sample_mesh(seed);
  double wgt = weight_.empty() ? 1.0 : weight_[elem_idx];
  return {u, wgt};
}

} // namespace openmc
