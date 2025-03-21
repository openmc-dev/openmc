#include "openmc/hex_mesh.h"
#include <algorithm> // for copy, equal, min, min_element
#include <cmath>     // for ceil, sqrt
#include <cstddef>   // for size_t
#include <gsl/gsl-lite.hpp>
#include <string>

#ifdef OPENMC_MPI
#include "mpi.h"
#endif

#include "xtensor/xbuilder.hpp"
#include "xtensor/xeval.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xnorm.hpp"
#include "xtensor/xsort.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include <fmt/core.h> // for fmt

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/memory.h"
#include "openmc/message_passing.h"
#include "openmc/openmp_interface.h"
#include "openmc/particle_data.h"
#include "openmc/plot.h"
#include "openmc/random_dist.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/volume_calc.h"
#include "openmc/xml_interface.h"

/*#ifdef LIBMESH*/
/*#include "libmesh/mesh_modification.h"*/
/*#include "libmesh/mesh_tools.h"*/
/*#include "libmesh/numeric_vector.h"*/
/*#endif*/

#ifdef DAGMC
#include "moab/FileOptions.hpp"
#endif

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

HexagonalMesh::HexagonalMesh(pugi::xml_node node)
  : PeriodicStructuredMesh {node}
{
  // Determine number of dimensions for mesh
  if (!check_for_node(node, "dimension")) {
    fatal_error("Must specify <dimension> on a regular mesh.");
  }

  xt::xtensor<int, 1> shape = get_node_xarray<int>(node, "dimension");
  int n = n_dimension_ = shape.size();
  if (n != 1 && n != 2) {
    fatal_error("Mesh must be one or two dimensions.");
  }
  std::copy(shape.begin(), shape.end(), shape_.begin());

  // Check that dimensions are all greater than zero
  if (xt::any(shape <= 0)) {
    fatal_error("All entries on the <dimension> element for a tally "
                "mesh must be positive.");
  }

  // Check for lower-left coordinates
  if (check_for_node(node, "lower_left")) {
    // Read mesh lower-left corner location
    lower_left_ = get_node_xarray<double>(node, "lower_left");
  } else {
    fatal_error("Must specify <lower_left> on a hexagonal mesh.");
  }

  // Make sure lower_left and dimension match
  if (shape.size() != lower_left_.size()) {
    fatal_error("Number of entries on <lower_left> must be the same "
                "as the number of entries on <dimension>.");
  }

  if (check_for_node(node, "width")) {
    // Make sure one of upper-right or width were specified
    if (check_for_node(node, "upper_right")) {
      fatal_error(
        "Cannot specify both <upper_right> and <width> on a hexgonal mesh.");
    }

    // Check to ensure width has same dimensions
    auto n = width_.size();
    if (n != lower_left_.size()) {
      fatal_error("Number of entries on <width> must be the same as "
                  "the number of entries on <lower_left>.");
    }

    // Check for negative widths
    if (xt::any(width_ < 0.0)) {
      fatal_error("Cannot have a negative <width> on a tally hexgonal mesh.");
    }

    // Set width and upper right coordinate
    upper_right_ = xt::eval(lower_left_ + shape * width_);

  } else if (check_for_node(node, "upper_right")) {
    upper_right_ = get_node_xarray<double>(node, "upper_right");

    // Check to ensure width has same dimensions
    auto n = upper_right_.size();
    if (n != lower_left_.size()) {
      fatal_error("Number of entries on <upper_right> must be the "
                  "same as the number of entries on <lower_left>.");
    }

    // Check that upper-right is above lower-left
    if (xt::any(upper_right_ < lower_left_)) {
      fatal_error("The <upper_right> coordinates must be greater than "
                  "the <lower_left> coordinates on a tally hexgonal mesh.");
    }

    // Set width
    width_ = xt::eval((upper_right_ - lower_left_) / shape);
  } else {
    fatal_error(
      "Must specify either <upper_right> or <width> on a hexagonal mesh.");
  }

  // n.b. must check that the number of hexes is odd
  //   or allow a parameter crown/ring
  int max_a = (shape[0] - 1) / 2;
  if (max_a == 0)
    hex_count_ = 1;
  else
    hex_count_ = 1 + 3 * (max_a) * (max_a - 1);

  // width[1] is the height of the full mesh block, width[0] is the width of
  // the full mesh block at its widest
  element_volume_ = width_[1] * width_[0] * width_[0] * sqrt(3);

  // Set material volumes
  volume_frac_ = 1.0 / hex_count_;

  // size of hex is defined as the radius of the circumscribed circle
  size_ = (width_[0] / shape[0]) / sqrt(3.0);

  // set the plane normals or 3 principal directions in the hex mash
  init_plane_normals();

  // scale basis vectors of hexagonal mesh
  scale_basis_vectors(size_);
}

const std::string HexagonalMesh::mesh_type = "hexagonal";

int HexagonalMesh::scale_basis_vectors(double s)
{
  // scale basis vectors of hexagonal mesh
  r_ = r_ * s;
  q_ = q_ * s;
  s_ = s_ * s;
  return 0;
}

int HexagonalMesh::init_plane_normals()
{
  n0_ = 0.5 * (r_ + q_ * 0.5);
  n1_ = 0.5 * (-s_ - r_ * 0.5);
  n2_ = 0.5 * (q_ + s_ * 0.5);
  return 0;
}

std::string HexagonalMesh::get_mesh_type() const
{
  return mesh_type;
}

int HexagonalMesh::hex_distance(
  const HexMeshIndex& ijkl0, const HexMeshIndex& ijkl1) const
{
  // return the integer lateral hex-distance between two hexes (ignores z)
  return std::max(
    std::max(std::abs(ijkl0[0] - ijkl1[0]), std::abs(ijkl0[1] - ijkl1[1])),
    std::abs(ijkl0[2] - ijkl1[2]));
}

int HexagonalMesh::hex_radius(const HexMeshIndex& ijkl) const
{
  // return the integer hex-radius of a hex. (ignores z)
  return std::max(
    std::max(std::abs(ijkl[0]), std::abs(ijkl[1])), std::abs(ijkl[2]));
}

int32_t HexagonalMesh::get_bin_from_hexindices(const HexMeshIndex& ijkl) const
{
  // get linear index from the HexMeshIndex
  int32_t r_0 = shape_[0];
  int32_t azim = 6 * (r_0 - 1);
  return ijkl[3] * hex_count_ + (1 + 3 * r_0 * (r_0 - 1)) + azim;
}

HexagonalMesh::HexMeshIndex HexagonalMesh::rotate_hexindex(
  const HexMeshIndex& ijkl, int steps) const
{
  HexMeshIndex new_ijkl {ijkl};
  for (int i; i < steps; i++) {
    int r = -new_ijkl[2];
    int q = -new_ijkl[0];
    int s = -new_ijkl[1];
    new_ijkl[0] = r;
    new_ijkl[1] = q;
    new_ijkl[2] = s;
  }
  return new_ijkl;
}

HexagonalMesh::HexMeshIndex HexagonalMesh::get_hexindices_from_bin(
  const int32_t bin) const
{
  // This does not yet take more than 1 z-slice into account
  int ring = (int)floor((sqrt(12 * bin - 3) + 3) / 6);

  int diff = bin - (1 + 3 * (ring - 1) * ring);
  HexMeshIndex ijkl {ring, 0, -ring, 0};
  return rotate_hexindex(ijkl, diff);
}

int HexagonalMesh::get_index_in_direction(const Position& r, int i) const
{
  switch (i) {
  case 0:
    return std::ceil((0.5 * r.x - 1.0 / (2 * sqrt(3)) * r.y) / this->size_);
  case 1:
    return std::ceil((1.0 / sqrt(3) * r.y) / this->size_);
  case 2:
    return std::ceil(
      ((0.5 * r.x - 1.0 / (2 * sqrt(3)) * r.y) - (1.0 / sqrt(3) * r.y)) /
      this->size_);
  case 3:
    // z
    return std::ceil((r.z - lower_left_[i]) / width_[i]);
  }
  return -1;
}

Position HexagonalMesh::get_position_from_hexindex(HexMeshIndex ijkl) const
{
  // return the cartesian position of center of a hexagon indexed by ijkl
  Position r;
  r.x = ijkl[0] * r_[0] * width_[0] + ijkl[1] * q_[0] * width_[0];
  r.y = ijkl[0] * r_[1] * width_[0] + ijkl[1] * q_[1] * width_[0];
  r.z = ijkl[3] * width_[1];

  return r;
}

HexagonalMesh::HexMeshIndex HexagonalMesh::get_hexindices(
  Position r, bool& in_mesh) const
{
  // return index of mesh element in hexes
  local_coords(r);

  HexMeshIndex idx;
  idx[0] = get_index_in_direction(r, 0);
  idx[1] = get_index_in_direction(r, 1);
  idx[2] = -idx[0] - idx[1];
  idx[3] = get_index_in_direction(r, 3);

  // check if either index is out of bounds
  in_mesh = in_hexmesh(idx);
  return idx;
}

bool HexagonalMesh::in_hexmesh(HexMeshIndex& ijkl) const
{
  for (auto it = ijkl.begin(); it < ijkl.end(); it++) {
    auto elem = *it;
    if (abs(elem) > shape_[0])
      return false;
  }
  return true;
}

Position HexagonalMesh::sample_element(
  const HexMeshIndex& ijkl, uint64_t* seed) const
{
  // return random position inside mesh element
  double r_r = uniform_distribution(0, 1, seed);
  double q_r = uniform_distribution(0, 1, seed);

  double x =
    this->r_[0] * (r_r + (ijkl[0] - 1)) + this->q_[0] * (q_r + (ijkl[1] - 1));
  double y =
    this->r_[1] * (r_r + (ijkl[0] - 1)) + this->q_[1] * (q_r + (ijkl[1] - 1));

  double z = uniform_distribution(-width_[1] / 2.0, width_[1] / 2.0, seed);

  return origin_ + Position(x, y, z);
}

double HexagonalMesh::find_r_crossing(
  const Position& r, const Direction& u, double l) const
{
  // solve r.x^2 + r.y^2 == r0^2
  // x^2 + 2*s*u*x + s^2*u^2 + s^2*v^2+2*s*v*y + y^2 -r0^2 = 0
  // s^2 * (u^2 + v^2) + 2*s*(u*x+v*y) + x^2+y^2-r0^2 = 0

  const double r0 = r_encl_;
  if (r0 == 0.0)
    return INFTY;

  const double denominator = u.x * u.x + u.y * u.y;

  // Direction of flight is in z-direction. Will never intersect r.
  if (std::abs(denominator) < FP_PRECISION)
    return INFTY;

  // inverse of dominator to help the compiler to speed things up
  const double inv_denominator = 1.0 / denominator;

  const double p = (u.x * r.x + u.y * r.y) * inv_denominator;
  double c = r.x * r.x + r.y * r.y - r0 * r0;
  double D = p * p - c * inv_denominator;

  if (D < 0.0)
    return INFTY;

  D = std::sqrt(D);

  // the solution -p - D is always smaller as -p + D : Check this one first
  if (std::abs(c) <= RADIAL_MESH_TOL)
    return INFTY;

  if (-p - D > l)
    return -p - D;
  if (-p + D > l)
    return -p + D;

  return INFTY;
}

template<class T>
void HexagonalMesh::raytrace_mesh(
  Position r0, Position r1, const Direction& u, T tally) const
{
  // TODO: when c++-17 is available, use "if constexpr ()" to compile-time
  // enable/disable tally calls for now, T template type needs to provide both
  // surface and track methods, which might be empty. modern optimizing
  // compilers will (hopefully) eliminate the complete code (including
  // calculation of parameters) but for the future: be explicit

  // very similar to to the structured mesh raytrace_mesh but we cannot rely on
  // simply recomputing only a single distance. Also the distance to the outside
  // of the mesh is somewhat more complicated

  // Compute the length of the entire track.
  double total_distance = (r1 - r0).norm();
  if (total_distance == 0.0 && settings::solver_type != SolverType::RANDOM_RAY)
    return;

  const int n = n_dimension_ * 2;

  // Flag if position is inside the mesh
  bool in_mesh;

  // Position is r = r0 + u * traveled_distance, start at r0
  double traveled_distance {0.0};

  // Calculate index of current cell. Offset the position a tiny bit in
  // direction of flight
  HexMeshIndex ijkl = get_indices(r0 + TINY_BIT * u, in_mesh);

  // if track is very short, assume that it is completely inside one cell.
  // Only the current cell will score and no surfaces
  if (total_distance < 2 * TINY_BIT) {
    if (in_mesh) {
      tally.track(ijkl, 1.0);
    }
    return;
  }

  // translate start and end positions,
  // this needs to come after the get_indices call because it does its own
  // translation
  local_coords(r0);
  local_coords(r1);

  // Calculate initial distances to next surfaces in all three dimensions
  std::array<MeshDistance, 8> distances;
  for (int k = 0; k < n; ++k) {
    distances[k] = distance_to_grid_boundary(ijkl, k, r0, u, 0.0);
  }

  // Loop until r = r1 is eventually reached
  while (true) {

    if (in_mesh) {

      // find surface with minimal distance to current position
      const auto k = std::min_element(distances.begin(), distances.end()) -
                     distances.begin();

      // Tally track length delta since last step
      tally.track(ijkl,
        (std::min(distances[k].distance, total_distance) - traveled_distance) /
          total_distance);

      // update position and leave, if we have reached end position
      traveled_distance = distances[k].distance;
      if (traveled_distance >= total_distance)
        return;

      // If we have not reached r1, we have hit a surface. Tally outward current
      tally.surface(ijkl, k, distances[k].max_surface, false);

      // Update cell and calculate distance to next surface in k-direction.
      ijkl = distances[k].next_index;
      // now the index has been updated recompute the distances - unfortunately
      // we have to do more than one again (as opposed to for cartesian mesh)

      if (k < 6) {
        for (int j = 0; j < 6; ++j) {
          distances[j] =
            distance_to_grid_boundary(ijkl, j, r0, u, traveled_distance);
        }
      } else {
        distances[6] =
          distance_to_grid_boundary(ijkl, 6, r0, u, traveled_distance);
      }
      // Check if we have left the interior of the mesh
      // Do this by getting new index
      in_mesh = true;
      for (auto i = 0; i < ijkl.size(); i++) {
        if (ijkl[i] < 1 || ijkl[i] > shape_[i]) {
          in_mesh = false;
          break;
        }
      }

      // If we are still inside the mesh, tally inward current for the next cell
      if (in_mesh)
        tally.surface(ijkl, k, !distances[k].max_surface, true);

    } else { // not inside mesh
      // we do this by the following algorithm:
      // 1. find a cylinder that completely circumscribes the hex-mesh region
      // 2. find the largest distance to the cylinder or to the z-plane
      // 3. add this to distance travelled
      //
      // If cylinder is used then do a single step minimal boundary find - this
      // would be enough to move us in-mesh.

      // For all directions outside the mesh, find the distance that we need to
      // travel to reach the next surface. Use the largest distance, as only
      // this will cross all outer surfaces.
      double dist_to_enclosing_cyl = find_r_crossing(r0, u, traveled_distance);
      for (int k = 0; k < n - 1; ++k) {
        distances[k] =
          distance_to_grid_boundary(ijkl, k, r0, u, dist_to_enclosing_cyl);
      }
      // we now need to compare the minimum of the 6 hex surface distances with
      // the z-surface distance and pick the longest of the two.
      int k_hex_min {0};
      int k_max {0};
      for (int k = 0; k < n - 1; ++k) {
        if (distances[k_hex_min].distance > distances[k].distance) {
          k_hex_min = k;
        }
      }
      if (distances[6].distance > distances[k_hex_min].distance) {
        k_max = 6;
      } else {
        k_max = k_hex_min;
      }
      traveled_distance = distances[k_max].distance;

      // If r1 is not inside the mesh, exit here
      if (traveled_distance >= total_distance)
        return;

      // Calculate the new cell index and update all distances to next surfaces.
      ijkl = get_indices(r0 + (traveled_distance + TINY_BIT) * u, in_mesh);
      for (int k = 0; k < n; ++k) {
        distances[k] =
          distance_to_grid_boundary(ijkl, k, r0, u, traveled_distance);
      }

      // If inside the mesh, Tally inward current
      if (in_mesh)
        tally.surface(ijkl, k_max, !distances[k_max].max_surface, true);
    }
  }
}

HexagonalMesh::HexMeshDistance HexagonalMesh::distance_to_hex_boundary(
  const HexMeshIndex& ijkl, int i, const Position& r0, const Direction& u,
  double l) const
{
  // Compute the distance to the element boundary of index i \in {0, ..., 6}
  // i==6 means z

  Position r = r0 - origin_;
  // Given the hex index - we now find the distance from r0 to the 0:q, 1:r, 2:s
  // plane(s) successively, given that the hex-center is given by the index
  // triplet and hence also the planes.
  Position rh = get_position_from_hexindex(ijkl);
  // local position relative to hex center
  Position rloc = r0 + l * u - rh;
  HexMeshDistance d;
  d.next_index = ijkl;

  if (i < 6) {

    switch (i) {
    case 0:
      d.distance = (this->r_ - rloc).dot(this->n0_) / u.dot(this->n0_);
      d.next_index[0]++;
      d.next_index[2]--;
      break;
    case 1:
      d.distance = (-this->s_ - rloc).dot(this->n1_) / u.dot(this->n1_);
      d.next_index[1]++;
      d.next_index[2]--;
      break;
    case 2:
      d.distance = (this->q_ - rloc).dot(this->n2_) / u.dot(this->n2_);
      d.next_index[0]--;
      d.next_index[1]++;
      break;
    case 3:
      d.distance = (-this->r_ - rloc).dot(this->n0_) / u.dot(this->n0_);
      d.next_index[0]--;
      d.next_index[2]++;
      break;
    case 4:
      d.distance = (this->s_ - rloc).dot(this->n1_) / u.dot(this->n1_);
      d.next_index[1]--;
      d.next_index[2]++;
      break;
    case 5:
      d.distance = (-this->q_ - rloc).dot(this->n2_) / u.dot(this->n2_);
      d.next_index[0]++;
      d.next_index[1]--;
      break;
    };
  } else {
    // Z-planes z-index has index 3 (nr 4) in ijkl.
    d.max_surface = (u[2] > 0);
    if (d.max_surface) {
      d.next_index[3]++;
      d.distance = ((lower_left_[2] + ijkl[3] * width_[2]) - r0[2]) / u[2];
    } else {
      d.next_index[3]--;
      d.distance =
        ((lower_left_[2] + (ijkl[3] - 1) * width_[2]) - r0[2]) / u[2];
    }
  }
  return d;
}

std::pair<vector<double>, vector<double>> HexagonalMesh::plot(
  Position plot_ll, Position plot_ur) const
{
  fatal_error("Plot of hexagonal Mesh not implemented");

  // Figure out which axes lie in the plane of the plot.
  array<vector<double>, 2> axis_lines;
  return {axis_lines[0], axis_lines[1]};
}

void HexagonalMesh::to_hdf5_inner(hid_t group) const
{
  hid_t mesh_group = create_group(group, "mesh " + std::to_string(id_));

  close_group(mesh_group);
}

double HexagonalMesh::volume(const HexMeshIndex& ijkl) const
{
  double zdiff = (upper_right_[2] - lower_left_[2]) / shape_[2];
  return 6 * sqrt(3.0) * 0.25 * size_ * size_ * zdiff;
}

} // namespace openmc
