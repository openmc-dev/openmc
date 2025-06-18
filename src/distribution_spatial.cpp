#include "openmc/distribution_spatial.h"

#include "openmc/error.h"
#include "openmc/mesh.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// SpatialDistribution implementation
//==============================================================================

unique_ptr<SpatialDistribution> SpatialDistribution::create(pugi::xml_node node)
{
  // Check for type of spatial distribution and read
  std::string type;
  if (check_for_node(node, "type"))
    type = get_node_value(node, "type", true, true);
  if (type == "cartesian") {
    return UPtrSpace {new CartesianIndependent(node)};
  } else if (type == "cylindrical") {
    return UPtrSpace {new CylindricalIndependent(node)};
  } else if (type == "spherical") {
    return UPtrSpace {new SphericalIndependent(node)};
  } else if (type == "mesh") {
    return UPtrSpace {new MeshSpatial(node)};
  } else if (type == "cloud") {
    return UPtrSpace {new PointCloud(node)};
  } else if (type == "box") {
    return UPtrSpace {new SpatialBox(node)};
  } else if (type == "fission") {
    return UPtrSpace {new SpatialBox(node, true)};
  } else if (type == "point") {
    return UPtrSpace {new SpatialPoint(node)};
  } else {
    fatal_error(fmt::format(
      "Invalid spatial distribution for external source: {}", type));
  }
}

//==============================================================================
// CartesianIndependent implementation
//==============================================================================

CartesianIndependent::CartesianIndependent(pugi::xml_node node)
{
  // Read distribution for x coordinate
  if (check_for_node(node, "x")) {
    pugi::xml_node node_dist = node.child("x");
    x_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at x=0
    double x[] {0.0};
    double p[] {1.0};
    x_ = UPtrDist {new Discrete {x, p, 1}};
  }

  // Read distribution for y coordinate
  if (check_for_node(node, "y")) {
    pugi::xml_node node_dist = node.child("y");
    y_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at y=0
    double x[] {0.0};
    double p[] {1.0};
    y_ = UPtrDist {new Discrete {x, p, 1}};
  }

  // Read distribution for z coordinate
  if (check_for_node(node, "z")) {
    pugi::xml_node node_dist = node.child("z");
    z_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at z=0
    double x[] {0.0};
    double p[] {1.0};
    z_ = UPtrDist {new Discrete {x, p, 1}};
  }
}

std::pair<Position, double> CartesianIndependent::sample(uint64_t* seed) const
{
  auto [x_val, x_wgt] = x_->sample(seed);
  auto [y_val, y_wgt] = y_->sample(seed);
  auto [z_val, z_wgt] = z_->sample(seed);
  Position xi {x_val, y_val, z_val};
  return { xi, x_wgt * y_wgt * z_wgt };
}

//==============================================================================
// CylindricalIndependent implementation
//==============================================================================

CylindricalIndependent::CylindricalIndependent(pugi::xml_node node)
{
  // Read distribution for r-coordinate
  if (check_for_node(node, "r")) {
    pugi::xml_node node_dist = node.child("r");
    r_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at r=0
    double x[] {0.0};
    double p[] {1.0};
    r_ = make_unique<Discrete>(x, p, 1);
  }

  // Read distribution for phi-coordinate
  if (check_for_node(node, "phi")) {
    pugi::xml_node node_dist = node.child("phi");
    phi_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at phi=0
    double x[] {0.0};
    double p[] {1.0};
    phi_ = make_unique<Discrete>(x, p, 1);
  }

  // Read distribution for z-coordinate
  if (check_for_node(node, "z")) {
    pugi::xml_node node_dist = node.child("z");
    z_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at z=0
    double x[] {0.0};
    double p[] {1.0};
    z_ = make_unique<Discrete>(x, p, 1);
  }

  // Read cylinder center coordinates
  if (check_for_node(node, "origin")) {
    auto origin = get_node_array<double>(node, "origin");
    if (origin.size() == 3) {
      origin_ = origin;
    } else {
      fatal_error(
        "Origin for cylindrical source distribution must be length 3");
    }
  } else {
    // If no coordinates were specified, default to (0, 0, 0)
    origin_ = {0.0, 0.0, 0.0};
  }
}

std::pair<Position, double> CylindricalIndependent::sample(uint64_t* seed) const
{
  auto [r, r_wgt] = r_->sample(seed);
  auto [phi, phi_wgt] = phi_->sample(seed);
  auto [z, z_wgt] = z_->sample(seed);
  double x = r * cos(phi) + origin_.x;
  double y = r * sin(phi) + origin_.y;
  z += origin_.z;
  Position xi {x, y, z};
  return { xi, r_wgt * phi_wgt * z_wgt};
}

//==============================================================================
// SphericalIndependent implementation
//==============================================================================

SphericalIndependent::SphericalIndependent(pugi::xml_node node)
{
  // Read distribution for r-coordinate
  if (check_for_node(node, "r")) {
    pugi::xml_node node_dist = node.child("r");
    r_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at r=0
    double x[] {0.0};
    double p[] {1.0};
    r_ = make_unique<Discrete>(x, p, 1);
  }

  // Read distribution for cos_theta-coordinate
  if (check_for_node(node, "cos_theta")) {
    pugi::xml_node node_dist = node.child("cos_theta");
    cos_theta_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at
    // cos_theta=0
    double x[] {0.0};
    double p[] {1.0};
    cos_theta_ = make_unique<Discrete>(x, p, 1);
  }

  // Read distribution for phi-coordinate
  if (check_for_node(node, "phi")) {
    pugi::xml_node node_dist = node.child("phi");
    phi_ = distribution_from_xml(node_dist);
  } else {
    // If no distribution was specified, default to a single point at phi=0
    double x[] {0.0};
    double p[] {1.0};
    phi_ = make_unique<Discrete>(x, p, 1);
  }

  // Read sphere center coordinates
  if (check_for_node(node, "origin")) {
    auto origin = get_node_array<double>(node, "origin");
    if (origin.size() == 3) {
      origin_ = origin;
    } else {
      fatal_error("Origin for spherical source distribution must be length 3");
    }
  } else {
    // If no coordinates were specified, default to (0, 0, 0)
    origin_ = {0.0, 0.0, 0.0};
  }
}

std::pair<Position, double> SphericalIndependent::sample(uint64_t* seed) const
{
  auto [r, r_wgt] = r_->sample(seed);
  auto [cos_theta, cos_theta_wgt] = cos_theta_->sample(seed);
  auto [phi, phi_wgt] = phi_->sample(seed);
  // sin(theta) by sin**2 + cos**2 = 1
  double x = r * std::sqrt(1 - cos_theta * cos_theta) * cos(phi) + origin_.x;
  double y = r * std::sqrt(1 - cos_theta * cos_theta) * sin(phi) + origin_.y;
  double z = r * cos_theta + origin_.z;
  Position xi {x, y, z};
  return { xi, r_wgt * cos_theta_wgt * phi_wgt };
}

//==============================================================================
// MeshSpatial implementation
//==============================================================================

MeshSpatial::MeshSpatial(pugi::xml_node node)
{

  if (get_node_value(node, "type", true, true) != "mesh") {
    fatal_error(fmt::format(
      "Incorrect spatial type '{}' for a MeshSpatial distribution"));
  }

  // No in-tet distributions implemented, could include distributions for the
  // barycentric coords Read in unstructured mesh from mesh_id value
  int32_t mesh_id = std::stoi(get_node_value(node, "mesh_id"));
  // Get pointer to spatial distribution
  mesh_idx_ = model::mesh_map.at(mesh_id);

  const auto mesh_ptr = model::meshes.at(mesh_idx_).get();

  check_element_types();

  size_t n_bins = this->n_sources();
  std::vector<double> strengths(n_bins, 1.0);

  // Create cdfs for sampling for an element over a mesh
  // Volume scheme is weighted by the volume of each tet
  // File scheme is weighted by an array given in the xml file
  if (check_for_node(node, "strengths")) {
    strengths = get_node_array<double>(node, "strengths");
    if (strengths.size() != n_bins) {
      fatal_error(
        fmt::format("Number of entries in the source strengths array {} does "
                    "not match the number of entities in mesh {} ({}).",
          strengths.size(), mesh_id, n_bins));
    }
  }

  if (get_node_value_bool(node, "volume_normalized")) {
    for (int i = 0; i < n_bins; i++) {
      strengths[i] *= this->mesh()->volume(i);
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

      if (get_node_value_bool(node, "volume_normalized")) {
        for (int i = 0; i < n_bins; i++) {
        bias_strengths[i] *= this->mesh()->volume(i);
        }
      }

      span<const double> b{bias_strengths};
      elem_idx_dist_.apply_bias(b);
    } else {
      fatal_error(
        fmt::format(
          "Bias node for mesh {} found without strengths array.", mesh_id));
    }
  }
}

MeshSpatial::MeshSpatial(int32_t mesh_idx, span<const double> strengths)
  : mesh_idx_(mesh_idx)
{
  check_element_types();
  elem_idx_dist_.assign(strengths);
}

void MeshSpatial::check_element_types() const
{
  const auto umesh_ptr = dynamic_cast<const UnstructuredMesh*>(this->mesh());
  if (umesh_ptr) {
    // ensure that the unstructured mesh contains only linear tets
    for (int bin = 0; bin < umesh_ptr->n_bins(); bin++) {
      if (umesh_ptr->element_type(bin) != ElementType::LINEAR_TET) {
        fatal_error(
          "Mesh specified for source must contain only linear tetrahedra.");
      }
    }
  }
}

int32_t MeshSpatial::sample_element_index(uint64_t* seed) const
{
  return elem_idx_dist_.sample(seed);
}

std::pair<int32_t, Position> MeshSpatial::sample_mesh(uint64_t* seed) const
{
  // Sample the CDF defined in initialization above
  int32_t elem_idx = this->sample_element_index(seed);
  return {elem_idx, mesh()->sample_element(elem_idx, seed)};
}

std::pair<Position, double> MeshSpatial::sample(uint64_t* seed) const
{
  auto [elem_idx, u] = this->sample_mesh(seed);
  return { u, elem_idx_dist_.weight()[elem_idx] };
}

//==============================================================================
// PointCloud implementation
//==============================================================================

PointCloud::PointCloud(pugi::xml_node node)
{
  if (check_for_node(node, "coords")) {
    point_cloud_ = get_node_position_array(node, "coords");
  } else {
    fatal_error("No coordinates were provided for the PointCloud "
                "spatial distribution");
  }

  std::vector<double> strengths;

  if (check_for_node(node, "strengths"))
    strengths = get_node_array<double>(node, "strengths");
  else
    strengths.resize(point_cloud_.size(), 1.0);

  if (strengths.size() != point_cloud_.size()) {
    fatal_error(
      fmt::format("Number of entries for the strengths array {} does "
                  "not match the number of spatial points provided {}.",
        strengths.size(), point_cloud_.size()));
  }

  point_idx_dist_.assign(strengths);

    if (check_for_node(node, "bias")) {
    pugi::xml_node bias_node = node.child("bias");

    if (check_for_node(bias_node, "strengths")) {
      std::vector<double> bias_strengths(point_cloud_.size(), 1.0);
      bias_strengths = get_node_array<double>(node, "strengths");

      if (bias_strengths.size() != point_cloud_.size()) {
        fatal_error(
          fmt::format("Number of entries in the bias strengths array {} does "
                      "not match the number of spatial points provided {}.",
            bias_strengths.size(), point_cloud_.size()));
      }

      span<const double> b{bias_strengths};
      point_idx_dist_.apply_bias(b);
    } else {
      fatal_error(
        fmt::format("Bias node for PointCloud found without strengths array."));
    }
  }
}

PointCloud::PointCloud(
  std::vector<Position> point_cloud, span<const double> strengths)
{
  point_cloud_.assign(point_cloud.begin(), point_cloud.end());
  point_idx_dist_.assign(strengths);
}

std::pair<Position, double> PointCloud::sample(uint64_t* seed) const
{
  int32_t index = point_idx_dist_.sample(seed);
  return { point_cloud_[index], point_idx_dist_.weight()[index] };
}

//==============================================================================
// SpatialBox implementation
//==============================================================================

SpatialBox::SpatialBox(pugi::xml_node node, bool fission)
  : only_fissionable_ {fission}
{
  // Read lower-right/upper-left coordinates
  auto params = get_node_array<double>(node, "parameters");
  if (params.size() != 6)
    openmc::fatal_error("Box/fission spatial source must have six "
                        "parameters specified.");

  lower_left_ = Position {params[0], params[1], params[2]};
  upper_right_ = Position {params[3], params[4], params[5]};
}

std::pair<Position, double> SpatialBox::sample(uint64_t* seed) const
{
  Position xi {prn(seed), prn(seed), prn(seed)};
  // Biasing not implemented--use CartesianIndependent instead
  return { lower_left_ + xi * (upper_right_ - lower_left_), 1.0 };
}

//==============================================================================
// SpatialPoint implementation
//==============================================================================

SpatialPoint::SpatialPoint(pugi::xml_node node)
{
  // Read location of point source
  auto params = get_node_array<double>(node, "parameters");
  if (params.size() != 3)
    openmc::fatal_error("Point spatial source must have three "
                        "parameters specified.");

  // Set position
  r_ = Position {params.data()};
}

std::pair<Position, double> SpatialPoint::sample(uint64_t* seed) const
{
  return { r_, 1.0 };
}

} // namespace openmc
