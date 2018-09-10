#include "openmc/mesh.h"

#include <algorithm> // for copy, min
#include <cstddef> // for size_t
#include <cmath>  // for ceil
#include <string>

#ifdef OPENMC_MPI
#include "mpi.h"
#endif
#include "xtensor/xbuilder.hpp"
#include "xtensor/xeval.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xsort.hpp"
#include "xtensor/xtensor.hpp"

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/message_passing.h"
#include "openmc/search.h"
#include "openmc/tallies/tally_filter.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

std::vector<std::unique_ptr<RegularMesh>> meshes;

std::unordered_map<int32_t, int32_t> mesh_map;

//==============================================================================
// RegularMesh implementation
//==============================================================================

RegularMesh::RegularMesh(pugi::xml_node node)
{
  // Copy mesh id
  if (check_for_node(node, "id")) {
    id_ = std::stoi(get_node_value(node, "id"));

    // Check to make sure 'id' hasn't been used
    if (mesh_map.find(id_) != mesh_map.end()) {
      fatal_error("Two or more meshes use the same unique ID: " +
        std::to_string(id_));
    }
  }

  // Read mesh type
  if (check_for_node(node, "type")) {
    auto temp = get_node_value(node, "type", true, true);
    if (temp == "regular") {
      // TODO: move elsewhere
    } else {
      fatal_error("Invalid mesh type: " + temp);
    }
  }

  // Determine number of dimensions for mesh
  if (check_for_node(node, "dimension")) {
    shape_ = get_node_xarray<int>(node, "dimension");
    int n = n_dimension_ = shape_.size();
    if (n != 1 && n != 2 && n != 3) {
      fatal_error("Mesh must be one, two, or three dimensions.");
    }

    // Check that dimensions are all greater than zero
    if (xt::any(shape_ <= 0)) {
      fatal_error("All entries on the <dimension> element for a tally "
        "mesh must be positive.");
    }
  }

  // Check for lower-left coordinates
  if (check_for_node(node, "lower_left")) {
    // Read mesh lower-left corner location
    lower_left_ = get_node_xarray<double>(node, "lower_left");
  } else {
    fatal_error("Must specify <lower_left> on a mesh.");
  }

  if (check_for_node(node, "width")) {
    // Make sure both upper-right or width were specified
    if (check_for_node(node, "upper_right")) {
      fatal_error("Cannot specify both <upper_right> and <width> on a mesh.");
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
      fatal_error("Cannot have a negative <width> on a tally mesh.");
    }

    // Set width and upper right coordinate
    upper_right_ = xt::eval(lower_left_ + shape_ * width_);

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
        "the <lower_left> coordinates on a tally mesh.");
    }

    // Set width and upper right coordinate
    width_ = xt::eval((upper_right_ - lower_left_) / shape_);
  } else {
    fatal_error("Must specify either <upper_right> and <width> on a mesh.");
  }

  if (shape_.dimension() > 0) {
    if (shape_.size() != lower_left_.size()) {
      fatal_error("Number of entries on <lower_left> must be the same "
        "as the number of entries on <dimension>.");
    }

    // Set volume fraction
    volume_frac_ = 1.0/xt::prod(shape_)();
  }
}

int RegularMesh::get_bin(Position r) const
{
  // Loop over the dimensions of the mesh
  for (int i = 0; i < n_dimension_; ++i) {
    // Check for cases where particle is outside of mesh
    if (r[i] < lower_left_[i]) {
      return -1;
    } else if (r[i] > upper_right_[i]) {
      return -1;
    }
  }

  // Determine indices
  int ijk[n_dimension_];
  bool in_mesh;
  get_indices(r, ijk, &in_mesh);
  if (!in_mesh) return -1;

  // Convert indices to bin
  return get_bin_from_indices(ijk);
}

int RegularMesh::get_bin_from_indices(const int* ijk) const
{
  if (n_dimension_ == 1) {
    return ijk[0];
  } else if (n_dimension_ == 2) {
    return (ijk[1] - 1)*shape_[0] + ijk[0];
  } else if (n_dimension_ == 3) {
    return ((ijk[2] - 1)*shape_[1] + (ijk[1] - 1))*shape_[0] + ijk[0];
  }
}

void RegularMesh::get_indices(Position r, int* ijk, bool* in_mesh) const
{
  // Find particle in mesh
  *in_mesh = true;
  for (int i = 0; i < n_dimension_; ++i) {
    ijk[i] = std::ceil((r[i] - lower_left_[i]) / width_[i]);

    // Check if indices are within bounds
    if (ijk[i] < 1 || ijk[i] > shape_[i]) *in_mesh = false;
  }
}

void RegularMesh::get_indices_from_bin(int bin, int* ijk) const
{
  if (n_dimension_ == 1) {
    ijk[0] = bin;
  } else if (n_dimension_ == 2) {
    ijk[0] = (bin - 1) % shape_[0] + 1;
    ijk[1] = (bin - 1) / shape_[0] + 1;
  } else if (n_dimension_ == 3) {
    ijk[0] = (bin - 1) % shape_[0] + 1;
    ijk[1] = ((bin - 1) % (shape_[0] * shape_[1])) / shape_[0] + 1;
    ijk[2] = (bin - 1) / (shape_[0] * shape_[1]) + 1;
  }
}

bool RegularMesh::intersects(Position r0, Position r1) const
{
  switch(n_dimension_) {
    case 1:
      return intersects_1d(r0, r1);
    case 2:
      return intersects_2d(r0, r1);
    case 3:
      return intersects_3d(r0, r1);
  }
}

bool RegularMesh::intersects_1d(Position r0, Position r1) const
{
  // Copy coordinates of mesh lower_left and upper_right
  double left = lower_left_[0];
  double right = upper_right_[0];

  // Check if line intersects either left or right surface
  if (r0.x < left) {
    return r1.x > left;
  } else if (r0.x < right) {
    return r1.x < left || r1.x > right;
  } else {
    return r1.x < right;
  }
}

bool RegularMesh::intersects_2d(Position r0, Position r1) const
{
  // Copy coordinates of starting point
  double x0 = r0.x;
  double y0 = r0.y;

  // Copy coordinates of ending point
  double x1 = r1.x;
  double y1 = r1.y;

  // Copy coordinates of mesh lower_left
  double xm0 = lower_left_[0];
  double ym0 = lower_left_[1];

  // Copy coordinates of mesh upper_right
  double xm1 = upper_right_[0];
  double ym1 = upper_right_[1];

  // Check if line intersects left surface -- calculate the intersection point y
  if ((x0 < xm0 && x1 > xm0) || (x0 > xm0 && x1 < xm0)) {
    double yi = y0 + (xm0 - x0) * (y1 - y0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1) {
      return true;
    }
  }

  // Check if line intersects back surface -- calculate the intersection point
  // x
  if ((y0 < ym0 && y1 > ym0) || (y0 > ym0 && y1 < ym0)) {
    double xi = x0 + (ym0 - y0) * (x1 - x0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1) {
      return true;
    }
  }

  // Check if line intersects right surface -- calculate the intersection
  // point y
  if ((x0 < xm1 && x1 > xm1) || (x0 > xm1 && x1 < xm1)) {
    double yi = y0 + (xm1 - x0) * (y1 - y0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1) {
      return true;
    }
  }

  // Check if line intersects front surface -- calculate the intersection point
  // x
  if ((y0 < ym1 && y1 > ym1) || (y0 > ym1 && y1 < ym1)) {
    double xi = x0 + (ym1 - y0) * (x1 - x0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1) {
      return true;
    }
  }
  return false;
}

bool RegularMesh::intersects_3d(Position r0, Position r1) const
{
  // Copy coordinates of starting point
  double x0 = r0.x;
  double y0 = r0.y;
  double z0 = r0.z;

  // Copy coordinates of ending point
  double x1 = r1.x;
  double y1 = r1.y;
  double z1 = r1.z;

  // Copy coordinates of mesh lower_left
  double xm0 = lower_left_[0];
  double ym0 = lower_left_[1];
  double zm0 = lower_left_[2];

  // Copy coordinates of mesh upper_right
  double xm1 = upper_right_[0];
  double ym1 = upper_right_[1];
  double zm1 = upper_right_[2];

  // Check if line intersects left surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm0 && x1 > xm0) || (x0 > xm0 && x1 < xm0)) {
    double yi = y0 + (xm0 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm0 - x0) * (z1 - z0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1 && zi >= zm0 && zi < zm1) {
      return true;
    }
  }

  // Check if line intersects back surface -- calculate the intersection point
  // (x,z)
  if ((y0 < ym0 && y1 > ym0) || (y0 > ym0 && y1 < ym0)) {
    double xi = x0 + (ym0 - y0) * (x1 - x0) / (y1 - y0);
    double zi = z0 + (ym0 - y0) * (z1 - z0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1 && zi >= zm0 && zi < zm1) {
      return true;
    }
  }

  // Check if line intersects bottom surface -- calculate the intersection
  // point (x,y)
  if ((z0 < zm0 && z1 > zm0) || (z0 > zm0 && z1 < zm0)) {
    double xi = x0 + (zm0 - z0) * (x1 - x0) / (z1 - z0);
    double yi = y0 + (zm0 - z0) * (y1 - y0) / (z1 - z0);
    if (xi >= xm0 && xi < xm1 && yi >= ym0 && yi < ym1) {
      return true;
    }
  }

  // Check if line intersects right surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm1 && x1 > xm1) || (x0 > xm1 && x1 < xm1)) {
    double yi = y0 + (xm1 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm1 - x0) * (z1 - z0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1 && zi >= zm0 && zi < zm1) {
      return true;
    }
  }

  // Check if line intersects front surface -- calculate the intersection point
  // (x,z)
  if ((y0 < ym1 && y1 > ym1) || (y0 > ym1 && y1 < ym1)) {
    double xi = x0 + (ym1 - y0) * (x1 - x0) / (y1 - y0);
    double zi = z0 + (ym1 - y0) * (z1 - z0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1 && zi >= zm0 && zi < zm1) {
      return true;
    }
  }

  // Check if line intersects top surface -- calculate the intersection point
  // (x,y)
  if ((z0 < zm1 && z1 > zm1) || (z0 > zm1 && z1 < zm1)) {
    double xi = x0 + (zm1 - z0) * (x1 - x0) / (z1 - z0);
    double yi = y0 + (zm1 - z0) * (y1 - y0) / (z1 - z0);
    if (xi >= xm0 && xi < xm1 && yi >= ym0 && yi < ym1) {
      return true;
    }
  }
  return false;
}

void RegularMesh::bins_crossed(const Particle* p, std::vector<int>& bins,
                               std::vector<double>& lengths) const
{
  constexpr int MAX_SEARCH_ITER = 100;

  // ========================================================================
  // Determine if the track intersects the tally mesh.

  // Copy the starting and ending coordinates of the particle.  Offset these
  // just a bit for the purposes of determining if there was an intersection
  // in case the mesh surfaces coincide with lattice/geometric surfaces which
  // might produce finite-precision errors.
  Position last_r {p->last_xyz};
  Position r {p->coord[0].xyz};
  Direction u {p->coord[0].uvw};

  Position r0 = last_r + TINY_BIT*u;
  Position r1 = r - TINY_BIT*u;

  // Determine indices for starting and ending location.
  int n = n_dimension_;
  xt::xtensor<int, 1> ijk0 = xt::empty<int>({n});
  bool start_in_mesh;
  get_indices(r0, ijk0.data(), &start_in_mesh);
  xt::xtensor<int, 1> ijk1 = xt::empty<int>({n});
  bool end_in_mesh;
  get_indices(r1, ijk1.data(), &end_in_mesh);

  // Check if the track intersects any part of the mesh.
  if (!start_in_mesh && !end_in_mesh) {
    if (!intersects(r0, r1)) return;
  }

  // ========================================================================
  // Figure out which mesh cell to tally.

  // Copy the un-modified coordinates the particle direction.
  r0 = last_r;
  r1 = r;

  // Compute the length of the entire track.
  double total_distance = (r1 - r0).norm();

  // We are looking for the first valid mesh bin.  Check to see if the
  // particle starts inside the mesh.
  if (!start_in_mesh) {
    xt::xtensor<double, 1> d = xt::zeros<double>({n});

    // The particle does not start in the mesh.  Note that we nudged the
    // start and end coordinates by a TINY_BIT each so we will have
    // difficulty resolving tracks that are less than 2*TINY_BIT in length.
    // If the track is that short, it is also insignificant so we can
    // safely ignore it in the tallies.
    if (total_distance < 2*TINY_BIT) return;

    // The particle does not start in the mesh so keep iterating the ijk0
    // indices to cross the nearest mesh surface until we've found a valid
    // bin.  MAX_SEARCH_ITER prevents an infinite loop.
    int search_iter = 0;
    int j;
    while (xt::any(ijk0 < 1) || xt::any(ijk0 > shape_)) {
      if (search_iter == MAX_SEARCH_ITER) {
        warning("Failed to find a mesh intersection on a tally mesh filter.");
        return;
      }

      for (j = 0; j < n; ++j) {
        if (std::fabs(u[j]) < FP_PRECISION) {
          d(j) = INFTY;
        } else if (u[j] > 0.0) {
          double xyz_cross = lower_left_[j] + ijk0(j) * width_[j];
          d(j) = (xyz_cross - r0[j]) / u[j];
        } else {
          double xyz_cross = lower_left_[j] + (ijk0(j) - 1) * width_[j];
          d(j) = (xyz_cross - r0[j]) / u[j];
        }
      }

      j = xt::argmin(d)(0);
      if (u[j] > 0.0) {
        ++ijk0(j);
      } else {
        --ijk0(j);
      }

      ++search_iter;
    }

    // Advance position
    r0 += d(j) * u;
  }

  while (true) {
    // ========================================================================
    // Compute the length of the track segment in the each mesh cell and return

    double distance;
    int j;
    if (ijk0 == ijk1) {
      // The track ends in this cell.  Use the particle end location rather
      // than the mesh surface.
      distance = (r1 - r0).norm();
    } else {
      // The track exits this cell.  Determine the distance to the closest mesh
      // surface.
      xt::xtensor<double, 1> d = xt::zeros<double>({n});
      for (int j = 0; j < n; ++j) {
        if (std::fabs(u[j]) < FP_PRECISION) {
          d(j) = INFTY;
        } else if (u[j] > 0) {
          double xyz_cross = lower_left_[j] + ijk0(j) * width_[j];
          d(j) = (xyz_cross - r0[j]) / u[j];
        } else {
          double xyz_cross = lower_left_[j] + (ijk0(j) - 1) * width_[j];
          d(j) = (xyz_cross - r0[j]) / u[j];
        }
      }
      j = xt::argmin(d)(0);
      distance = d(j);
    }

    // Assign the next tally bin and the score.
    int bin = get_bin_from_indices(ijk0.data());
    bins.push_back(bin);
    lengths.push_back(distance / total_distance);

    // If the particle track ends in that bin, then we are done.
    if (ijk0 == ijk1) break;

    // Translate the starting coordintes by the distance to that face. This
    // should be the xyz that we computed the distance to in the last
    // iteration of the filter loop.
    r0 += distance * u;

    // Increment the indices into the next mesh cell.
    if (u[j] > 0.0) {
      ++ijk0(j);
    } else {
      --ijk0(j);
    }

    // If the next indices are invalid, then the track has left the mesh and
    // we are done.
    if (xt::any(ijk0 < 1) || xt::any(ijk0 > shape_)) break;
  }
}

void RegularMesh::surface_bins_crossed(const Particle* p, std::vector<int>& bins) const
{
  // ========================================================================
  // Determine if the track intersects the tally mesh.

  // Copy the starting and ending coordinates of the particle.
  Position r0 {p->last_xyz_current};
  Position r1 {p->coord[0].xyz};
  Direction u {p->coord[0].uvw};

  // Determine indices for starting and ending location.
  int n = n_dimension_;
  xt::xtensor<int, 1> ijk0 = xt::empty<int>({n});
  bool start_in_mesh;
  get_indices(r0, ijk0.data(), &start_in_mesh);
  xt::xtensor<int, 1> ijk1 = xt::empty<int>({n});
  bool end_in_mesh;
  get_indices(r1, ijk1.data(), &end_in_mesh);

  // Check if the track intersects any part of the mesh.
  if (!start_in_mesh && !end_in_mesh) {
    if (!intersects(r0, r1)) return;
  }

  // ========================================================================
  // Figure out which mesh cell to tally.

  // Calculate number of surface crossings
  int n_cross = xt::sum(xt::abs(ijk1 - ijk0))();
  if (n_cross == 0) return;

  // Bounding coordinates
  Position xyz_cross;
  for (int i = 0; i < n; ++i) {
    if (u[i] > 0.0) {
      xyz_cross[i] = lower_left_[i] + ijk0[i] * width_[i];
    } else {
      xyz_cross[i] = lower_left_[i] + (ijk0[i] - 1) * width_[i];
    }
  }

  for (int j = 0; j < n_cross; ++j) {
    // Set the distances to infinity
    Position d {INFTY, INFTY, INFTY};

    // Determine closest bounding surface. We need to treat
    // special case where the cosine of the angle is zero since this would
    // result in a divide-by-zero.
    double distance = INFTY;
    for (int i = 0; i < n; ++i) {
      if (u[i] == 0) {
        d[i] = INFINITY;
      } else {
        d[i] = (xyz_cross[i] - r0[i])/u[i];
      }
      distance = std::min(distance, d[i]);
    }

    // Loop over the dimensions
    for (int i = 0; i < n; ++i) {
      // Check whether distance is the shortest distance
      if (distance == d[i]) {

        // Check whether particle is moving in positive i direction
        if (u[i] > 0) {

          // Outward current on i max surface
          if (xt::all(ijk0 >= 1) && xt::all(ijk0 <= shape_)) {
            int i_surf = 4*i + 3;
            int i_mesh = get_bin_from_indices(ijk0.data());
            int i_bin = 4*n*(i_mesh - 1) + i_surf;

            bins.push_back(i_bin);
          }

          // Advance position
          ++ijk0[i];
          xyz_cross[i] += width_[i];

          // If the particle crossed the surface, tally the inward current on
          // i min surface
          if (xt::all(ijk0 >= 1) && xt::all(ijk0 <= shape_)) {
            int i_surf = 4*i + 2;
            int i_mesh = get_bin_from_indices(ijk0.data());
            int i_bin = 4*n*(i_mesh - 1) + i_surf;

            bins.push_back(i_bin);
          }

        } else {
          // The particle is moving in the negative i direction

          // Outward current on i min surface
          if (xt::all(ijk0 >= 1) && xt::all(ijk0 <= shape_) ){
            int i_surf = 4*i + 1;
            int i_mesh = get_bin_from_indices(ijk0.data());
            int i_bin = 4*n*(i_mesh - 1) + i_surf;

            bins.push_back(i_bin);
          }

          // Advance position
          --ijk0[i];
          xyz_cross[i] -= width_[i];

          // If the particle crossed the surface, tally the inward current on
          // i max surface
          if (xt::all(ijk0 >= 1) && xt::all(ijk0 <= shape_)) {
            int i_surf = 4*i + 4;
            int i_mesh = get_bin_from_indices(ijk0.data());
            int i_bin = 4*n*(i_mesh - 1) + i_surf;

            bins.push_back(i_bin);
          }
        }
      }
    }

    // Calculate new coordinates
    r0 += distance * u;
  }
}

void RegularMesh::to_hdf5(hid_t group) const
{
  hid_t mesh_group = create_group(group, "mesh " + std::to_string(id_));

  write_dataset(mesh_group, "type", "regular");
  write_dataset(mesh_group, "dimension", shape_);
  write_dataset(mesh_group, "lower_left", lower_left_);
  write_dataset(mesh_group, "upper_right", upper_right_);
  write_dataset(mesh_group, "width", width_);

  close_group(mesh_group);
}

xt::xarray<double> RegularMesh::count_sites(int64_t n, const Bank* bank,
  int n_energy, const double* energies, bool* outside) const
{
  // Determine shape of array for counts
  std::size_t m = xt::prod(shape_)();
  std::vector<std::size_t> shape;
  if (n_energy > 0) {
    shape = {m, static_cast<std::size_t>(n_energy - 1)};
  } else {
    shape = {m};
  }

  // Create array of zeros
  xt::xarray<double> cnt {shape, 0.0};
  bool outside_ = false;

  for (int64_t i = 0; i < n; ++i) {
    // determine scoring bin for entropy mesh
    // TODO: off-by-one
    int mesh_bin = get_bin({bank[i].xyz}) - 1;

    // if outside mesh, skip particle
    if (mesh_bin < 0) {
      outside_ = true;
      continue;
    }

    if (n_energy > 0) {
      double E = bank[i].E;
      if (E >= energies[0] && E <= energies[n_energy - 1]) {
        // determine energy bin
        int e_bin = lower_bound_index(energies, energies + n_energy, E);

        // Add to appropriate bin
        cnt(mesh_bin, e_bin) += bank[i].wgt;
      }
    } else {
      // Add to appropriate bin
      cnt(mesh_bin) += bank[i].wgt;
    }
  }

  // Create copy of count data
  int total = cnt.size();
  double* cnt_reduced = new double[total];

#ifdef OPENMC_MPI
  // collect values from all processors
  MPI_Reduce(cnt.data(), cnt_reduced, total, MPI_DOUBLE, MPI_SUM, 0,
    mpi::intracomm);

  // Check if there were sites outside the mesh for any processor
  if (outside) {
    MPI_Reduce(&outside_, outside, 1, MPI_C_BOOL, MPI_LOR, 0, mpi::intracomm);
  }
#else
  std::copy(cnt.data(), cnt.data() + total, cnt_reduced);
  if (outside) *outside = outside_;
#endif

  // Adapt reduced values in array back into an xarray
  auto arr = xt::adapt(cnt_reduced, total, xt::acquire_ownership(), shape);
  xt::xarray<double> counts = arr;

  return counts;
}

//==============================================================================
// C API functions
//==============================================================================

//! Extend the meshes array by n elements
extern "C" int
openmc_extend_meshes(int32_t n, int32_t* index_start, int32_t* index_end)
{
  if (index_start) *index_start = meshes.size();
  for (int i = 0; i < n; ++i) {
    meshes.emplace_back(new RegularMesh{});
  }
  if (index_end) *index_end = meshes.size() - 1;

  return 0;
}

//! Return the index in the meshes array of a mesh with a given ID
extern "C" int
openmc_get_mesh_index(int32_t id, int32_t* index)
{
  auto pair = mesh_map.find(id);
  if (pair == mesh_map.end()) {
    set_errmsg("No mesh exists with ID=" + std::to_string(id) + ".");
    return OPENMC_E_INVALID_ID;
  }
  *index = pair->second;
  return 0;
}

// Return the ID of a mesh
extern "C" int
openmc_mesh_get_id(int32_t index, int32_t* id)
{
  if (index < 0 || index >= meshes.size()) {
    set_errmsg("Index in meshes array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  *id = meshes[index]->id_;
  return 0;
}

//! Set the ID of a mesh
extern "C" int
openmc_mesh_set_id(int32_t index, int32_t id)
{
  if (index < 0 || index >= meshes.size()) {
    set_errmsg("Index in meshes array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  meshes[index]->id_ = id;
  mesh_map[id] = index;
  return 0;
}

//! Get the dimension of a mesh
extern "C" int
openmc_mesh_get_dimension(int32_t index, int** dims, int* n)
{
  if (index < 0 || index >= meshes.size()) {
    set_errmsg("Index in meshes array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  *dims = meshes[index]->shape_.data();
  *n = meshes[index]->n_dimension_;
  return 0;
}

//! Set the dimension of a mesh
extern "C" int
openmc_mesh_set_dimension(int32_t index, int n, const int* dims)
{
  if (index < 0 || index >= meshes.size()) {
    set_errmsg("Index in meshes array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  // Copy dimension
  std::vector<std::size_t> shape = {static_cast<std::size_t>(n)};
  auto& m = meshes[index];
  m->shape_ = xt::adapt(dims, n, xt::no_ownership(), shape);
  m->n_dimension_ = m->shape_.size();

  return 0;
}

//! Get the mesh parameters
extern "C" int
openmc_mesh_get_params(int32_t index, double** ll, double** ur, double** width, int* n)
{
  if (index < 0 || index >= meshes.size()) {
    set_errmsg("Index in meshes array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  auto& m = meshes[index];
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

//! Set the mesh parameters
extern "C" int
openmc_mesh_set_params(int32_t index, int n, const double* ll, const double* ur,
                       const double* width)
{
  if (index < 0 || index >= meshes.size()) {
    set_errmsg("Index in meshes array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  auto& m = meshes[index];
  std::vector<std::size_t> shape = {static_cast<std::size_t>(n)};
  if (ll && ur) {
    m->lower_left_ = xt::adapt(ll, n, xt::no_ownership(), shape);
    m->upper_right_ = xt::adapt(ur, n, xt::no_ownership(), shape);
    m->width_ = (m->upper_right_ - m->lower_left_) / m->shape_;
  } else if (ll && width) {
    m->lower_left_ = xt::adapt(ll, n, xt::no_ownership(), shape);
    m->width_ = xt::adapt(width, n, xt::no_ownership(), shape);
    m->upper_right_ = m->lower_left_ + m->shape_ * m->width_;
  } else if (ur && width) {
    m->upper_right_ = xt::adapt(ur, n, xt::no_ownership(), shape);
    m->width_ = xt::adapt(width, n, xt::no_ownership(), shape);
    m->lower_left_ = m->upper_right_ - m->shape_ * m->width_;
  } else {
    set_errmsg("At least two parameters must be specified.");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  return 0;
}

//==============================================================================
// Non-member functions
//==============================================================================

void read_meshes(pugi::xml_node* root)
{
  for (auto node : root->children("mesh")) {
    // Read mesh and add to vector
    meshes.emplace_back(new RegularMesh{node});

    // Map ID to position in vector
    mesh_map[meshes.back()->id_] = meshes.size() - 1;
  }
}

void meshes_to_hdf5(hid_t group)
{
  // Write number of meshes
  hid_t meshes_group = create_group(group, "meshes");
  int32_t n_meshes = meshes.size();
  write_attribute(meshes_group, "n_meshes", n_meshes);

  if (n_meshes > 0) {
    // Write IDs of meshes
    std::vector<int> ids;
    for (const auto& m : meshes) {
      m->to_hdf5(meshes_group);
      ids.push_back(m->id_);
    }
    write_attribute(meshes_group, "ids", ids);
  }

  close_group(meshes_group);
}

//==============================================================================
// Fortran compatibility
//==============================================================================

extern "C" {
  int n_meshes() { return meshes.size(); }

  RegularMesh* mesh_ptr(int i) { return meshes.at(i).get(); }

  int32_t mesh_id(RegularMesh* m) { return m->id_; }

  double mesh_volume_frac(RegularMesh* m) { return m->volume_frac_; }

  int mesh_n_dimension(RegularMesh* m) { return m->n_dimension_; }

  int mesh_dimension(RegularMesh* m, int i) { return m->shape_(i - 1); }

  double mesh_lower_left(RegularMesh* m, int i) { return m->lower_left_(i - 1); }

  double mesh_upper_right(RegularMesh* m, int i) { return m->upper_right_(i - 1); }

  double mesh_width(RegularMesh* m, int i) { return m->width_(i - 1); }

  int mesh_get_bin(RegularMesh* m, const double* xyz)
  {
    return m->get_bin({xyz});
  }

  int mesh_get_bin_from_indices(RegularMesh* m, const int* ijk)
  {
    return m->get_bin_from_indices(ijk);
  }

  void mesh_get_indices(RegularMesh* m, const double* xyz, int* ijk, bool* in_mesh)
  {
    m->get_indices({xyz}, ijk, in_mesh);
  }

  void mesh_get_indices_from_bin(RegularMesh* m, int bin, int* ijk)
  {
    m->get_indices_from_bin(bin, ijk);
  }

  void mesh_bins_crossed(RegularMesh* m, const Particle* p,
    TallyFilterMatch* match)
  {
    match->bins.clear();
    match->weights.clear();
    m->bins_crossed(p, match->bins, match->weights);
  }

  void mesh_surface_bins_crossed(RegularMesh* m, const Particle* p,
    TallyFilterMatch* match)
  {
    match->bins.clear();
    match->weights.clear();
    m->surface_bins_crossed(p, match->bins);
    for (auto b : match->bins) {
      match->weights.push_back(1.0);
    }
  }

  void free_memory_mesh()
  {
    meshes.clear();
    mesh_map.clear();
  }
}


} // namespace openmc
