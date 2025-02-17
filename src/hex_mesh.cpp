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
#include "xtensor/xsort.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xnorm.hpp"
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

#ifdef LIBMESH
const bool LIBMESH_ENABLED = true;
#else
const bool LIBMESH_ENABLED = false;
#endif

namespace model {

std::unordered_map<int32_t, int32_t> mesh_map;
vector<unique_ptr<Mesh>> meshes;

} // namespace model

HexagonalMesh::HexagonalMesh(pugi::xml_node node) : StructuredMesh {node}
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
      fatal_error("Cannot specify both <upper_right> and <width> on a hexgonal mesh.");
    }

    width_ = get_node_xarray<double>(node, "width");

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
    fatal_error("Must specify either <upper_right> or <width> on a hexagonal mesh.");
  }

  // n.b. must check that the number of hexes is odd
  //   or allow a parameter crown/ring
  int max_a_ = (shape[0] - 1) / 2;
  if (max_a == 0)
    hex_count_ = 1;
  else
    hex_count_ = 1 + 3 * (max_a) * (max_a-1);

  // width[1] is the height of the full mesh block, width[0] is the width of
  // the full mesh block at its widest
  element_volume_ = width_[1] * width_[0] * width_[0] * sqrt(3);

  // Set material volumes
  volume_frac_ = 1.0 / hex_count_;

  // size of hex is defined as the radius of the circumscribed circle
  size_ = (width_[0]/shape[0])/sqrt(3.0);

  // set the plane normals or 3 principal directions in the hex mash
  init_plane_normals();

  // scale basis vectors of hexagonal mesh
  scale_basis_vectors(size_);
}

const std::string Hegagonal::mesh_type = "hexagonal";

int HexagonalMesh::scale_basis_vectors(double s){
  // scale basis vectors of hexagonal mesh
  a_ = a_ * s;
  b_ = b_ * s;
  c_ = c_ * s;
  return 0;
}

int HaxgonalMesh::init_plane_normals(){
  n0_ = 0.5 * (a_ - c_);
  n1_ = 0.5 * (b_ - c_);
  n2_ = 0.5 * (-a_ + c_);
  return 0;
}

std::string HexgonalMesh::get_mesh_type() const
{
  return mesh_type;
}

int HexagonalMesh::hex_distance(const HexMesIndex& ijkl0, const HexMesIndex& ijkl1) const
{
  // return the integer lateral hex-distance between two hexes (ignores z)
  return std::max( std::max(std::abs(ijkl0[0]-ijkl1[0]), std::abs(ijkl0[1]-ijkl1[1])), std::abs(ijkl0[2]-ijkl1[2]) );
}

int HexagonalMesh::hex_radius(const HexMeshIndex &ijkl) const
{
  //return the integer hex-radius of a hex. (ignores z)
  return std::max( std::max(std::abs(ijkl[0]), std::abs(ijkl[1])), std::abs(ijkl[2]) );
}

int HexgonalMesh::get_bin_from_indices(const HexMeshIndex& ijkl) const
{
  //get linear index from the HexMeshIndex
  auto r_0 = hex_radius();
  auto azim = 6 * (r_0-1);
  return ijkl[3] * hex_count_ + (1 + 3*r0*(r0-1)) + azim;
}

int HexagonalMesh::get_index_in_direction(const Position& r, int i) const
{
  switch (i){
    case 0:
      return std::ceil( (0.5 * r.x - 1.0/(2*sqrt(3)) * r.y) / this->abs_a );
    case 1:
      return std::ceil( (1.0/sqrt(3) * r.y) / this->abs_a );
    case 2:
      return std::ceil( ( (0.5 * r.x - 1.0/(2*sqrt(3)) * r.y) - (1.0/sqrt(3) * r.y) )/ this->abs_a );
    case 3:
      // z
      return std::ceil((r - lower_left_[i]) / width_[i]);
  }
}

Position HexagonalMesh::get_positition_of_hexindex(HexMeshIndex ijkl) const
{
  //return the cartesian position of a hex indexed by abcz
  Position r;
  r.x = ijkl[0]*a_[0]*width_[0] + ijkl[1]*b_[0]*width_[0];
  r.y = ijkl[0]*a_[1]*width_[0] + ijkl[1]*b_[1]*width_[0];;
  r.z = ijkl[3]*width_[1];
}

HexagonalMesh::HexMeshIndex HexgonalMesh::get_indices(
  Position r, bool& in_mesh) const
{
  //return index of mesh element in hexes
  local_coords(r);

  HexMeshIndex idx;
  idx[0] = get_index_in_direction(r,0);
  idx[1] = get_index_in_direction(r,1);
  idx[2] = -idx[0] -idx[1];
  idx[3] = get_index_in_direction(r, 3);

  //check if either index is out of bounds
  in_mesh = in_bounds();
  return idx;
}

bool HexagonalMesh::in_bounds(HexMeshIndex& ijkl) const
{
  for (auto it=ijkl.begin(); it<ijkl.end(); it++){
    auto elem = *it;
    if (abs(elem)>max_a_) return false;
  }
  return true;
}


Position HexgonalMesh::sample_element(
  const HexMeshIndex& ijk, uint64_t* seed) const
{
  //return random position inside mesh element
  double a_r = uniform_distribution();
  double b_r = uniform_distribution();

  x = this->a_[0] * (a_r + (ijk[0]-1) ) + this->b_[0] * (b_r + (ijk[1]-1) );
  y = this->a_[1] * (a_r + (ijk[0]-1) ) + this->b_[1] * (b_r + (ijk[1]-1) );

  double z = uniform_distribution(-width_[1]/2.0, width_[1]/2.0, seed);

  return origin_ + Position(x, y, z);
}


double HexgonalMesh::find_r_crossing(
  const Position& r, const Direction& u, double l, int shell) const
{
  //finds the cylinder radius where r+u crosses a plane
  //maybe relevant to hexagonal since it wil tell us approximately which
  //hex-ring we're in
  if ((shell < 0) || (shell > shape_[0]))
    return INFTY;

  // solve r.x^2 + r.y^2 == r0^2
  // x^2 + 2*s*u*x + s^2*u^2 + s^2*v^2+2*s*v*y + y^2 -r0^2 = 0
  // s^2 * (u^2 + v^2) + 2*s*(u*x+v*y) + x^2+y^2-r0^2 = 0

  const double r0 = grid_[0][shell];
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

StructuredMesh::MeshDistance HexgonalMesh::find_z_crossing(
  const Position& r, const Direction& u, double l, int shell) const
{
  //finds plane-crossing in z - should be the same as for Cylindrical.

  MeshDistance d;
  d.next_index = shell;

  // Direction of flight is within xy-plane. Will never intersect z.
  if (std::abs(u.z) < FP_PRECISION)
    return d;

  d.max_surface = (u.z > 0.0);
  if (d.max_surface && (shell <= shape_[2])) {
    d.next_index += 1;
    d.distance = (grid_[2][shell] - r.z) / u.z;
  } else if (!d.max_surface && (shell > 0)) {
    d.next_index -= 1;
    d.distance = (grid_[2][shell - 1] - r.z) / u.z;
  }
  return d;
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
  // simply recomputing only a single distance. Also the distance to the outside of the mesh
  // is somewhat more complicated

  // Compute the length of the entire track.
  double total_distance = (r1 - r0).norm();
  if (total_distance == 0.0 && settings::solver_type != SolverType::RANDOM_RAY)
    return;

  const int n = n_dimension_;

  // Flag if position is inside the mesh
  bool in_mesh;

  // Position is r = r0 + u * traveled_distance, start at r0
  double traveled_distance {0.0};

  // Calculate index of current cell. Offset the position a tiny bit in
  // direction of flight
  MeshIndex ijk = get_indices(r0 + TINY_BIT * u, in_mesh);

  // if track is very short, assume that it is completely inside one cell.
  // Only the current cell will score and no surfaces
  if (total_distance < 2 * TINY_BIT) {
    if (in_mesh) {
      tally.track(ijk, 1.0);
    }
    return;
  }

  // translate start and end positions,
  // this needs to come after the get_indices call because it does its own
  // translation
  local_coords(r0);
  local_coords(r1);

  // Calculate initial distances to next surfaces in all three dimensions
  std::array<MeshDistance, 3> distances;
  for (int k = 0; k < n; ++k) {
    distances[k] = distance_to_grid_boundary(ijk, k, r0, u, 0.0);
  }

  // Loop until r = r1 is eventually reached
  while (true) {

    if (in_mesh) {

      // find surface with minimal distance to current position
      const auto k = std::min_element(distances.begin(), distances.end()) -
                     distances.begin();

      // Tally track length delta since last step
      tally.track(ijk,
        (std::min(distances[k].distance, total_distance) - traveled_distance) /
          total_distance);

      // update position and leave, if we have reached end position
      traveled_distance = distances[k].distance;
      if (traveled_distance >= total_distance)
        return;

      // If we have not reached r1, we have hit a surface. Tally outward current
      tally.surface(ijk, k, distances[k].max_surface, false);

      // Update cell and calculate distance to next surface in k-direction.
      // The two other directions are still valid!
      ijk[k] = distances[k].next_index;
      distances[k] =
        distance_to_grid_boundary(ijk, k, r0, u, traveled_distance);

      // Check if we have left the interior of the mesh
      in_mesh = ((ijk[k] >= 1) && (ijk[k] <= shape_[k]));

      // If we are still inside the mesh, tally inward current for the next cell
      if (in_mesh)
        tally.surface(ijk, k, !distances[k].max_surface, true);

    } else { // not inside mesh

      // For all directions outside the mesh, find the distance that we need to
      // travel to reach the next surface. Use the largest distance, as only
      // this will cross all outer surfaces.
      int k_max {0};
      for (int k = 0; k < n; ++k) {
        if ((ijk[k] < 1 || ijk[k] > shape_[k]) &&
            (distances[k].distance > traveled_distance)) {
          traveled_distance = distances[k].distance;
          k_max = k;
        }
      }

      // If r1 is not inside the mesh, exit here
      if (traveled_distance >= total_distance)
        return;

      // Calculate the new cell index and update all distances to next surfaces.
      ijk = get_indices(r0 + (traveled_distance + TINY_BIT) * u, in_mesh);
      for (int k = 0; k < n; ++k) {
        distances[k] =
          distance_to_grid_boundary(ijk, k, r0, u, traveled_distance);
      }

      // If inside the mesh, Tally inward current
      if (in_mesh)
        tally.surface(ijk, k_max, !distances[k_max].max_surface, true);
    }
  }
}


StructuredMesh::MeshDistance HexgonalMesh::distance_to_grid_boundary(
  const MeshIndex& ijk, int i, const Position& r0, const Direction& u,
  double l) const
{
  // Compute the distance to the element boundary of index i
  // Since MeshIndex is defined to have 3 and only 3 numbers (can we override?) we need to do a bit of gymnastics to
  // convert to a HexMeshIndex
  //
  Position r = r0 - origin_;

  if (i < 4) {
    // we chack the hexagonal planes
    // have to check them all since we may have a new situtation after a boundary crossing.
  } else {
    // check z-planes
  }
}

double HaxegonalMesh::distance_to_hex_boundary(HexMeshIndex &ijkl, int i){
  // return the distance to the hex boundary in direction i
  // i==0: a, i==1: b, i==2: c
  // This is not independent of the index triplet ijkl
  MeshDistance d;
  switch (i){
    case 0:
      d.next_index=i;
      // next index cannot be used as a single scalar - we need to recompute all distances
      // compute the distance to the SE/n0_ plane
    case 1:
  }
  //
}



int HexagonalMesh::set_grid()
{
  // not necessary to do this - since it is a regular mesh we can do it on the fly
}

MeshDistance find_hex_crossing(const Position& r0, const direction& u, double l, const HexMeshIndex& ijkl, int i) const
{
  // return the distance to the hex-plane in the desired direction
  // i is the index number in which direction we'll be looking


}


std::pair<vector<double>, vector<double>> HexagonalMesh::plot(
  Position plot_ll, Position plot_ur) const
{
  fatal_error("Plot of hexagonal Mesh not implemented");

  // Figure out which axes lie in the plane of the plot.
  array<vector<double>, 2> axis_lines;
  return {axis_lines[0], axis_lines[1]};
}

void HexgonalMesh::to_hdf5(hid_t group) const
{
  hid_t mesh_group = create_group(group, "mesh " + std::to_string(id_));

  write_dataset(mesh_group, "type", "hexgonal");
  write_dataset(mesh_group, "a_grid", grid_[0]);
  write_dataset(mesh_group, "b_grid", grid_[1]);
  write_dataset(mesh_group, "z_grid", grid_[2]);
  write_dataset(mesh_group, "origin", origin_);

  close_group(mesh_group);
}

double HexgonalMesh::volume(const MeshIndex& ijk) const
{
  double a_i = grid_[0][ijk[0] - 1];
  double a_o = grid_[0][ijk[0]];

  double b_i = grid_[1][ijk[1] - 1];
  double b_o = grid_[1][ijk[1]];

  double z_i = grid_[2][ijk[2] - 1];
  double z_o = grid_[2][ijk[2]];

  return 6 * sqrt(3.0) * (a_o-a_i)*(b_o-b_i)*(z_o - z_i);
}
