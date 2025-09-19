#include "openmc/hex_mesh.h"
#include <algorithm> // for copy, equal, min, min_element
#include <cmath>     // for ceil, sqrt
#include <cstddef>   // for size_t
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

  // Make sure shape has two numbers
  if (n == 1) {
    shape_[1] = 1;
    n = 2;
  }

  // Check that dimensions are all greater than zero
  if (xt::any(shape <= 0)) {
    fatal_error("All entries on the <dimension> element for a tally "
                "mesh must be positive.");
  }
  if (shape_[0] % 2 == 0) {
    fatal_error("First shape coordinate must be odd to avoid non-integer "
                "ring count.");
  }
  hex_radius_ = (shape_[0] - 1) / 2;

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
    width_ = get_node_xarray<double>(node, "width");
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
  if (hex_radius_ == 0)
    hex_count_ = 1;
  else
    hex_count_ = 1 + 3 * (hex_radius_ + 1) * hex_radius_;

  // width[1] is the height of the full mesh block, width[0] is the width of
  // the hexagon from flat end to flat end
  element_volume_ = width_[1] * width_[0] * width_[0] * sqrt(3) * 0.5;

  // Set material volumes
  volume_frac_ = 1.0 / hex_count_ / shape_[1];

  // size of hex is defined as the radius of the circumscribed circle
  size_ = width_[0] / sqrt(3.0);

  // radius of enclosing cylinder
  r_encl_ = (hex_radius_ - 0.5) * sqrt(3) * size_ + (1 - sqrt(3) * 0.5) * size_;

  // set the plane normals or 3 principal directions in the hex mesh
  init_plane_normals();

  // scale grid vectors of hexagonal mesh
  scale_grid_vectors(size_);
}

const std::string HexagonalMesh::mesh_type = "hexagonal";

double HexagonalMesh::volume(const StructuredMesh::MeshIndex& ijk) const
{
  return element_volume_;
}

int HexagonalMesh::scale_grid_vectors(double s)
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

  n0_ /= n0_.norm();
  n1_ /= n1_.norm();
  n2_ /= n2_.norm();

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
  int32_t r_0 = hex_radius(ijkl);
  if (r_0 == 0)
    return (ijkl[3] - 1) * hex_count_;
  int32_t start_of_ring = (1 + 3 * r_0 * (r_0 - 1));
  int32_t bin_no = (ijkl[3] - 1) * hex_count_ + (1 + 3 * r_0 * (r_0 - 1)) +
                   offset_in_ring(ijkl, r_0);
  return bin_no;
}

int32_t HexagonalMesh::offset_in_ring(const HexMeshIndex& ijkl, int32_t r) const
{
  // find the offset within a ring
  if (r == 0)
    return 0;
  HexMeshIndex i {ijkl};
  HexMeshIndex corner {r, 0, -r, 0};

  int segment {0};
  if (abs(i[2]) >= abs(i[1]) && abs(i[2]) >= abs(i[0])) {
    if (i[2] < 0)
      segment = 0;
    else
      segment = 3;
  } else if (abs(i[1]) >= abs(i[0])) {
    if (i[1] > 0)
      segment = 1;
    else
      segment = 4;
  } else {
    if (i[0] < 0)
      segment = 2;
    else
      segment = 5;
  }

  for (int k = 0; k < segment; k++)
    corner = rotate_hexindex(corner);

  int32_t hexd = hex_distance(corner, i);
  int ii = r * segment + hexd;
  return ii;
}

HexagonalMesh::HexMeshIndex HexagonalMesh::rotate_hexindex(
  const HexMeshIndex& ijkl) const
{
  HexMeshIndex new_ijkl {ijkl};
  new_ijkl[0] = -ijkl[1];
  new_ijkl[1] = -ijkl[2];
  new_ijkl[2] = -ijkl[0];
  return new_ijkl;
}

HexagonalMesh::HexMeshIndex HexagonalMesh::get_hexindices_from_bin(
  const int32_t bin) const
{
  std::array<HexMeshIndex, 6> directions;
  directions[0] = {-1, 1, 0, 0};
  directions[1] = {-1, 0, 1, 0};
  directions[2] = {0, -1, 1, 0};
  directions[3] = {1, -1, 0, 0};
  directions[4] = {1, 0, -1, 0};
  directions[5] = {0, 1, -1, 0};

  HexMeshIndex ijkl = {0, 0, 0, 1};
  ijkl[3] = (int)floor(bin / hex_count_) + 1;
  int spiral_index = bin % hex_count_;
  if (spiral_index == 0) {
    return ijkl;
  }
  int ring = (int)floor((sqrt(12 * spiral_index - 3) + 3) / 6);
  int start_of_ring = (1 + 3 * ring * (ring - 1));

  if (ring > 0) {
    int segment = (spiral_index - start_of_ring) / ring;

    ijkl[0] = ring;
    ijkl[2] = -ring;

    for (int k = 0; k < segment; k++)
      ijkl = rotate_hexindex(ijkl);

    for (int k = 0; k < spiral_index - start_of_ring - ring * segment; k++) {
      for (int l = 0; l < ijkl.size(); l++)
        ijkl[l] = ijkl[l] + directions[segment][l];
    }
  }
  return ijkl;
}

int HexagonalMesh::n_bins() const
{
  return hex_count_ * shape_[1];
}

int HexagonalMesh::n_surface_bins() const
{
  return 4 * 4 * n_bins();
}

xt::xtensor<int, 1> HexagonalMesh::get_x_shape() const
{
  // because method is const, shape_ is const as well and can't be adapted
  auto tmp_shape = shape_;
  return xt::adapt(tmp_shape, {2});
}

std::string HexagonalMesh::bin_label(int bin) const
{
  HexMeshIndex ijkl = get_hexindices_from_bin(bin);
  int hr = hex_radius(ijkl);
  int ofr = offset_in_ring(ijkl, hr);
  return fmt::format(
    "Mesh Index ({}, {}, {})", hr, offset_in_ring(ijkl, hr), ijkl[3]);
}

double HexagonalMesh::frac_hexindex_in_direction(const Position& r, int i) const
{
  switch (i) {
  case 0:
    return (2.0 / 3.0 * -r.y) / this->size_;
  case 1:
    return (sqrt(3.0) / 3.0 * r.x + 1.0 / 3.0 * r.y) / this->size_;
  case 2:
    return -(2.0 / 3.0 * -r.y) / this->size_ -
           (sqrt(3.0) / 3.0 * r.x + 1.0 / 3.0 * r.y) / this->size_;
  case 3:
    // z is idx 1 in width_ and lower_left_ / upper_right_
    return (r.z - lower_left_[1]) / width_[1];
  }
  return -1;
}

int HexagonalMesh::get_index_in_direction(double r, int i) const
{
  // dummy placeholder
  return 0;
}

Position HexagonalMesh::get_position_from_hexindex(HexMeshIndex ijkl) const
{
  // return the cartesian position of center of a hexagon indexed by ijkl
  // Note thet we have to use the plane normals as basis vectors, as opposed to
  // the grid vectors (r, q, s)
  Position r;
  r.x = ijkl[0] * n0_[0] * size_ * sqrt(3) + ijkl[1] * n1_[0] * size_ * sqrt(3);
  r.y = ijkl[0] * n0_[1] * size_ * sqrt(3) + ijkl[1] * n1_[1] * size_ * sqrt(3);
  r.z = (ijkl[3] - 0.5) * width_[1] + lower_left_[1];

  return r;
}

HexagonalMesh::HexMeshIndex HexagonalMesh::get_hexindices(
  Position r, bool& in_mesh) const
{
  // return index of mesh element in hexes
  local_coords(r);
  vector<double> frac_cds {0, 0, 0, 0};

  // r coordinate
  frac_cds[0] = frac_hexindex_in_direction(r, 0);
  // q coordinate
  frac_cds[1] = frac_hexindex_in_direction(r, 1);
  // s coordinate
  frac_cds[2] = frac_hexindex_in_direction(r, 2);
  // z-coordinate
  frac_cds[3] = frac_hexindex_in_direction(r, 3);

  HexMeshIndex idx = round_frac_hexindex(frac_cds);
  // check if either index is out of bounds
  in_mesh = in_hexmesh(idx);
  return idx;
}

HexagonalMesh::HexMeshIndex HexagonalMesh::round_frac_hexindex(
  vector<double> frac_ijkl) const
{
  std::vector<double> diff(4);
  HexMeshIndex ijkl {0, 0, 0, 1};

  for (int i = 0; i < frac_ijkl.size(); ++i) {
    diff[i] = (std::abs(std::round(frac_ijkl[i]) - frac_ijkl[i]));
    ijkl[i] = std::round(frac_ijkl[i]);
  }
  if (diff[0] > diff[1] && diff[0] > diff[2]) {
    ijkl[0] = -ijkl[1] - ijkl[2];
  } else if (diff[1] > diff[2]) {
    ijkl[1] = -ijkl[0] - ijkl[2];
  } else {
    ijkl[2] = -ijkl[0] - ijkl[1];
  }
  // z-coordinate should be treated differently
  ijkl[3] = std::ceil(frac_ijkl[3]);

  return ijkl;
}

bool HexagonalMesh::in_hexmesh(HexMeshIndex& ijkl) const
{
  for (auto it = ijkl.begin(); it != std::prev(ijkl.end()); it++) {
    auto elem = *it;
    if (abs(elem) > hex_radius_)
      return false;
  }
  if (ijkl[3] > shape_[1] || ijkl[3] <= 0)
    return false;
  return true;
}

Position HexagonalMesh::sample_element(
  const HexMeshIndex& ijkl, uint64_t* seed) const
{
  // return random position inside mesh element
  Position rh = get_position_from_hexindex(ijkl);

  Position xy = sample_hexagon(seed) * size_;
  Position z(
    0, 0, uniform_distribution(-width_[1] / 2.0, width_[1] / 2.0, seed));
  return origin_ + rh + xy + z;
}

Position HexagonalMesh::sample_hexagon(uint64_t* seed) const
{
  // return a position inside hexagon at (0,0,0)
  double x = uniform_distribution(-sqrt(3) * 0.5, sqrt(3) * 0.5, seed);
  double y = uniform_distribution(-0.5, 1.0, seed);

  if (y > 0.5)
    if (y - 0.5 > 0.5 - std::abs(x) / sqrt(3)) { // reflect the point
      if (x > 0)
        x -= sqrt(3) * 0.5;
      else
        x += sqrt(3) * 0.5;
      y -= 1.5;
    }
  return Position(x, y, 0);
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
  int iteration {0};
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
  HexMeshIndex ijkl = get_hexindices(r0 + TINY_BIT * u, in_mesh);

  // If outside mesh and not on the way towards it (i.e. missing the
  // circumscribed cylinder) exit early.
  double dist_to_enclosing_cyl = find_r_crossing(r0, u, traveled_distance);

  if (!in_mesh) {
    if (dist_to_enclosing_cyl < INFTY &&
        dist_to_enclosing_cyl < total_distance) {
      traveled_distance = dist_to_enclosing_cyl;
    } else
      return;
  }

  // if track is very short, assume that it is completely inside one cell.
  // Only the current cell will score and no surfaces
  if (total_distance < 2 * TINY_BIT) {
    if (in_mesh) {
      tally.track(ijkl, 1.0);
    }
    return;
  }

  // translate start and end positions,
  // this needs to come after the get_hexindices call because it does its own
  // translation
  local_coords(r0);
  local_coords(r1);

  // Calculate initial distances to next surfaces in all three dimensions
  std::array<HexMeshDistance, 4> distances;
  for (int k = 0; k < n; ++k) {
    distances[k] = distance_to_hex_boundary(ijkl, k, r0, u, traveled_distance);
  }

  // Loop until r = r1 is eventually reached
  while (true) {
    iteration++;
    // std::cout << iteration << std::endl;
    //  find surface with minimal distance to current position
    const auto k =
      std::min_element(distances.begin(), distances.end()) - distances.begin();

    if (in_mesh) {
      // Tally track length delta since last step
      tally.track(ijkl,
        (std::min(distances[k].distance, total_distance) - traveled_distance) /
          total_distance);
    }

    // update position and leave, if we have reached end position
    traveled_distance = distances[k].distance;
    if (traveled_distance >= total_distance)
      return;

    if (in_mesh) {
      // If we have not reached r1, we have hit a surface. Tally outward current
      tally.surface(ijkl, k, distances[k].max_surface, false);
    }

    // Update cell and calculate distance to next surface in k-direction.
    ijkl = distances[k].next_index;
    // now the index has been updated recompute the distances - unfortunately
    // we have to do more than one again (as opposed to for cartesian mesh)

    if (k < 3) {
      for (int j = 0; j < 4; ++j) {
        distances[j] =
          distance_to_hex_boundary(ijkl, j, r0, u, traveled_distance);
      }
    } else {
      distances[3] =
        distance_to_hex_boundary(ijkl, 3, r0, u, traveled_distance);
      for (int j = 0; j < 3; ++j) {
        distances[j].next_index[3] = ijkl[3];
      }
    }

    // Check if we have left the interior of the mesh
    // Do this by getting new index
    in_mesh = in_hexmesh(ijkl);

    // If we are still inside the mesh, tally inward current for the next cell
    if (in_mesh)
      tally.surface(ijkl, k, !distances[k].max_surface, true);
  }
}

HexagonalMesh::HexMeshDistance HexagonalMesh::distance_to_hex_boundary(
  const HexMeshIndex& ijkl, int i, const Position& r0, const Direction& u,
  double l) const
{
  // Compute the distance to the element boundary of index i \in {0, ..., 3}
  // i==3 means z

  Position r = r0 - origin_;
  // Given the hex index - we now find the distance from r0 to the 0:q, 1:r, 2:s
  // plane(s) successively, given that the hex-center is given by the hexindex
  // quadruplet and hence also the planes.
  Position rh = get_position_from_hexindex(ijkl);
  // local position relative to hex center
  Position rloc = r0 + l * u - rh;
  HexMeshDistance d;
  d.next_index = ijkl;

  double dh = 0;

  if (i < 3) {
    if (std::abs(u[0]) + std::abs(u[1]) < FP_PRECISION)
      return d;
    switch (i) {
    case 0:
      if (std::abs(u.dot(n0_)) < FP_PRECISION) {
        return d;
      }
      dh = rh.dot(n0_) / u.dot(n0_);
      d.max_surface = n0_.dot(u) > 0;
      if (abs(ijkl[0]) > hex_radius_ + 1) {
        return d;
      }
      if (d.max_surface) {
        d.distance = dh + (this->r_ - r0).dot(this->n0_) / u.dot(this->n0_);
        d.next_index[0]++;
        d.next_index[2]--;
      } else {
        d.distance = dh + (-this->r_ - r0).dot(this->n0_) / u.dot(this->n0_);
        d.next_index[0]--;
        d.next_index[2]++;
      }
      break;
    case 1:
      if (std::abs(u.dot(n1_)) < FP_PRECISION) {
        return d;
      }
      dh = rh.dot(n1_) / u.dot(n1_);
      d.max_surface = n1_.dot(u) > 0;
      if (abs(ijkl[1]) > hex_radius_ + 1) {
        return d;
      }
      if (d.max_surface) {
        d.distance = dh + (-this->s_ - r0).dot(this->n1_) / u.dot(this->n1_);
        d.next_index[1]++;
        d.next_index[2]--;
      } else {
        d.distance = dh + (this->s_ - r0).dot(this->n1_) / u.dot(this->n1_);
        d.next_index[1]--;
        d.next_index[2]++;
      }
      break;
    case 2:
      if (std::abs(u.dot(n2_)) < FP_PRECISION) {
        return d;
      }
      dh = rh.dot(n2_) / u.dot(n2_);
      d.max_surface = n2_.dot(u) > 0;
      if (abs(ijkl[2]) > hex_radius_ + 1) {
        return d;
      }
      if (d.max_surface) {
        d.distance = dh + (this->q_ - r0).dot(this->n2_) / u.dot(this->n2_);
        d.next_index[0]--;
        d.next_index[1]++;
      } else {
        d.distance = dh + (-this->q_ - r0).dot(this->n2_) / u.dot(this->n2_);
        d.next_index[0]++;
        d.next_index[1]--;
      }
      break;
    }
  } else {
    if (std::abs(u[2]) < FP_PRECISION) {
      return d;
    }
    // Z-planes z-index has index 3 (nr 4) in ijkl.
    d.max_surface = (u[2] > 0);
    if (d.max_surface) {
      d.next_index[3]++;
      d.distance = ((lower_left_[1] + ijkl[3] * width_[1]) - r0[2]) / u[2];
    } else {
      d.next_index[3]--;
      d.distance =
        ((lower_left_[1] + (ijkl[3] - 1) * width_[1]) - r0[2]) / u[2];
    }
  }
  return d;
}

StructuredMesh::MeshDistance HexagonalMesh::distance_to_grid_boundary(
  const MeshIndex& ijk, int i, const Position& r0, const Direction& u,
  double l) const
{
  MeshDistance d;
  d.distance = INFTY;
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

void HexagonalMesh::to_hdf5_inner(hid_t mesh_group) const
{
  write_dataset(mesh_group, "dimension", get_x_shape());
  write_dataset(mesh_group, "lower_left", lower_left_);
  write_dataset(mesh_group, "upper_right", upper_right_);
  write_dataset(mesh_group, "width", width_);
  // hex-mesh specifics
  write_dataset(mesh_group, "hex_count", hex_count_);
  write_dataset(mesh_group, "hex_radius", hex_radius_);
}

double HexagonalMesh::volume(const HexMeshIndex& ijkl) const
{
  double zdiff = (upper_right_[2] - lower_left_[2]) / shape_[2];
  return 6 * sqrt(3.0) * 0.25 * size_ * size_ * zdiff;
}

void HexagonalMesh::bins_crossed(Position r0, Position r1, const Direction& u,
  vector<int>& bins, vector<double>& lengths) const
{

  // Helper tally class.
  // stores a pointer to the mesh class and references to bins and lengths
  // parameters. Performs the actual tally through the track method.
  struct TrackAggregator {
    TrackAggregator(
      const HexagonalMesh* _mesh, vector<int>& _bins, vector<double>& _lengths)
      : mesh(_mesh), bins(_bins), lengths(_lengths)
    {}
    void surface(const HexMeshIndex& ijkl, int k, bool max, bool inward) const
    {}
    void track(const HexMeshIndex& ijkl, double l) const
    {
      bins.push_back(mesh->get_bin_from_hexindices(ijkl));
      lengths.push_back(l);
    }

    const HexagonalMesh* mesh;
    vector<int>& bins;
    vector<double>& lengths;
  };

  // Perform the mesh raytrace with the helper class.
  raytrace_mesh(r0, r1, u, TrackAggregator(this, bins, lengths));
}

void HexagonalMesh::surface_bins_crossed(
  Position r0, Position r1, const Direction& u, vector<int>& bins) const
{

  // Helper tally class.
  // stores a pointer to the mesh class and a reference to the bins parameter.
  // Performs the actual tally through the surface method.
  struct SurfaceAggregator {
    SurfaceAggregator(const HexagonalMesh* _mesh, vector<int>& _bins)
      : mesh(_mesh), bins(_bins)
    {}
    void surface(const HexMeshIndex& ijkl, int k, bool max, bool inward) const
    {
      int i_bin = 4 * 4 * mesh->get_bin_from_hexindices(ijkl) + 4 * k;
      if (max)
        i_bin += 2;
      if (inward)
        i_bin += 1;
      bins.push_back(i_bin);
    }
    void track(const HexMeshIndex& idx, double l) const {}

    const HexagonalMesh* mesh;
    vector<int>& bins;
  };

  // Perform the mesh raytrace with the helper class.
  raytrace_mesh(r0, r1, u, SurfaceAggregator(this, bins));
}

//==============================================================================
// Helper functions for the C API
//==============================================================================

int check_hex_mesh(int32_t index)
{
  if (index < 0 || index >= model::meshes.size()) {
    set_errmsg("Index in meshes array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  return 0;
}
// This is identical to the one in mesh.cpp
template<class T>
int check_mesh_type(int32_t index)
{
  if (int err = check_hex_mesh(index))
    return err;

  T* mesh = dynamic_cast<T*>(model::meshes[index].get());
  if (!mesh) {
    set_errmsg("This function is not valid for input mesh.");
    return OPENMC_E_INVALID_TYPE;
  }
  return 0;
}

// This is identical to the one in mesh.cpp
template<class T>
bool is_mesh_type(int32_t index)
{
  T* mesh = dynamic_cast<T*>(model::meshes[index].get());
  return mesh;
}

//! Get the dimension of a hexagonal mesh
extern "C" int openmc_hexagonal_mesh_get_dimension(
  int32_t index, int** dims, int* n)
{
  if (int err = check_mesh_type<HexagonalMesh>(index))
    return err;
  HexagonalMesh* mesh =
    dynamic_cast<HexagonalMesh*>(model::meshes[index].get());
  *dims = mesh->shape_.data();
  *n = mesh->n_dimension_;
  return 0;
}

//! Set the dimension of a hexagonal mesh
extern "C" int openmc_hexagonal_mesh_set_dimension(
  int32_t index, int n, const int* dims)
{
  if (int err = check_mesh_type<HexagonalMesh>(index))
    return err;
  HexagonalMesh* mesh =
    dynamic_cast<HexagonalMesh*>(model::meshes[index].get());

  // Copy dimension
  mesh->n_dimension_ = n;
  std::copy(dims, dims + n, mesh->shape_.begin());

  // TODO: incorporate this bit into method in HexagonalMesh that can be called
  // from here and from constructor
  mesh->hex_radius_ = (mesh->shape_[0] - 1) / 2;
  if (mesh->hex_radius_ == 0)
    mesh->hex_count_ = 1;
  else
    mesh->hex_count_ = 1 + 3 * (mesh->hex_radius_ + 1) * mesh->hex_radius_;

  return 0;
}
//! Get the hexagonal mesh parameters
extern "C" int openmc_hexagonal_mesh_get_params(
  int32_t index, double** ll, double** ur, double** width, int* n)
{
  if (int err = check_mesh_type<HexagonalMesh>(index))
    return err;
  HexagonalMesh* m = dynamic_cast<HexagonalMesh*>(model::meshes[index].get());

  if (m->lower_left_.dimension() == 0) {
    set_errmsg("Mesh parameters have not been set.");
    return OPENMC_E_ALLOCATE;
  }

  *ll = m->lower_left_.data();
  *ur = m->upper_right_.data();
  *width = m->width_.data();
  *n = m->n_dimension_;
  return 0;
}

//! Set the hexagonal mesh parameters
extern "C" int openmc_hexagonal_mesh_set_params(
  int32_t index, int n, const double* ll, const double* ur, const double* width)
{
  if (int err = check_mesh_type<HexagonalMesh>(index))
    return err;
  HexagonalMesh* m = dynamic_cast<HexagonalMesh*>(model::meshes[index].get());

  if (m->n_dimension_ == -1) {
    set_errmsg("Need to set mesh dimension before setting parameters.");
    return OPENMC_E_UNASSIGNED;
  }

  vector<std::size_t> shape = {static_cast<std::size_t>(n)};
  if (ll && ur) {
    m->lower_left_ = xt::adapt(ll, n, xt::no_ownership(), shape);
    m->upper_right_ = xt::adapt(ur, n, xt::no_ownership(), shape);
    m->width_ = (m->upper_right_ - m->lower_left_) / m->get_x_shape();
  } else if (ll && width) {
    m->lower_left_ = xt::adapt(ll, n, xt::no_ownership(), shape);
    m->width_ = xt::adapt(width, n, xt::no_ownership(), shape);
    m->upper_right_ = m->lower_left_ + m->get_x_shape() * m->width_;
  } else if (ur && width) {
    m->upper_right_ = xt::adapt(ur, n, xt::no_ownership(), shape);
    m->width_ = xt::adapt(width, n, xt::no_ownership(), shape);
    m->lower_left_ = m->upper_right_ - m->get_x_shape() * m->width_;
  } else {
    set_errmsg("At least two parameters must be specified.");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  // TODO: incorporate this into method in HexagonalMesh that can be called from
  // here and from constructor

  // Set material volumes etc.
  m->size_ = m->width_[0] / sqrt(3.0);
  m->init_plane_normals();
  m->scale_grid_vectors(m->size_);

  m->volume_frac_ = 1.0 / m->shape_[1] / m->hex_count_;
  m->element_volume_ =
    m->width_[1] * m->width_[0] * m->width_[0] * sqrt(3) * 0.5;

  return 0;
}

} // namespace openmc
