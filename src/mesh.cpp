#include "openmc/mesh.h"
#include <algorithm> // for copy, equal, min, min_element
#include <cmath>     // for ceil
#include <cstddef>   // for size_t
#include <gsl/gsl-lite.hpp>
#include <string>

#ifdef OPENMC_MPI
#include "mpi.h"
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#include "xtensor/xbuilder.hpp"
#include "xtensor/xeval.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xsort.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include <fmt/core.h> // for fmt

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/memory.h"
#include "openmc/message_passing.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/xml_interface.h"

#ifdef LIBMESH
#include "libmesh/mesh_tools.h"
#include "libmesh/numeric_vector.h"
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

#ifdef LIBMESH
namespace settings {
unique_ptr<libMesh::LibMeshInit> libmesh_init;
const libMesh::Parallel::Communicator* libmesh_comm {nullptr};
} // namespace settings
#endif

//==============================================================================
// Helper functions
//==============================================================================

//! Update an intersection point if the given candidate is closer.
//
//! The first 6 arguments are coordinates for the starting point of a particle
//! and its intersection with a mesh surface.  If the distance between these
//! two points is shorter than the given `min_distance`, then the `r` argument
//! will be updated to match the intersection point, and `min_distance` will
//! also be updated.

inline bool check_intersection_point(double x1, double x0, double y1, double y0,
  double z1, double z0, Position& r, double& min_distance)
{
  double dist =
    std::pow(x1 - x0, 2) + std::pow(y1 - y0, 2) + std::pow(z1 - z0, 2);
  if (dist < min_distance) {
    r.x = x1;
    r.y = y1;
    r.z = z1;
    min_distance = dist;
    return true;
  }
  return false;
}

//==============================================================================
// Mesh implementation
//==============================================================================

Mesh::Mesh(pugi::xml_node node)
{
  // Copy mesh id
  if (check_for_node(node, "id")) {
    id_ = std::stoi(get_node_value(node, "id"));

    // Check to make sure 'id' hasn't been used
    if (model::mesh_map.find(id_) != model::mesh_map.end()) {
      fatal_error(
        "Two or more meshes use the same unique ID: " + std::to_string(id_));
    }
  }
}

void Mesh::set_id(int32_t id)
{
  Expects(id >= 0 || id == C_NONE);

  // Clear entry in mesh map in case one was already assigned
  if (id_ != C_NONE) {
    model::mesh_map.erase(id_);
    id_ = C_NONE;
  }

  // Ensure no other mesh has the same ID
  if (model::mesh_map.find(id) != model::mesh_map.end()) {
    throw std::runtime_error {
      fmt::format("Two meshes have the same ID: {}", id)};
  }

  // If no ID is specified, auto-assign the next ID in the sequence
  if (id == C_NONE) {
    id = 0;
    for (const auto& m : model::meshes) {
      id = std::max(id, m->id_);
    }
    ++id;
  }

  // Update ID and entry in the mesh map
  id_ = id;
  model::mesh_map[id] = model::meshes.size() - 1;
}

//==============================================================================
// Structured Mesh implementation
//==============================================================================

std::string StructuredMesh::bin_label(int bin) const
{
  vector<int> ijk(n_dimension_);
  get_indices_from_bin(bin, ijk.data());

  if (n_dimension_ > 2) {
    return fmt::format("Mesh Index ({}, {}, {})", ijk[0], ijk[1], ijk[2]);
  } else if (n_dimension_ > 1) {
    return fmt::format("Mesh Index ({}, {})", ijk[0], ijk[1]);
  } else {
    return fmt::format("Mesh Index ({})", ijk[0]);
  }
}

//==============================================================================
// Unstructured Mesh implementation
//==============================================================================

UnstructuredMesh::UnstructuredMesh(pugi::xml_node node) : Mesh(node)
{

  // check the mesh type
  if (check_for_node(node, "type")) {
    auto temp = get_node_value(node, "type", true, true);
    if (temp != "unstructured") {
      fatal_error(fmt::format("Invalid mesh type: {}", temp));
    }
  }

  // get the filename of the unstructured mesh to load
  if (check_for_node(node, "filename")) {
    filename_ = get_node_value(node, "filename");
    if (!file_exists(filename_)) {
      fatal_error("Mesh file '" + filename_ + "' does not exist!");
    }
  } else {
    fatal_error(fmt::format(
      "No filename supplied for unstructured mesh with ID: {}", id_));
  }

  // check if mesh tally data should be written with
  // statepoint files
  if (check_for_node(node, "output")) {
    output_ = get_node_value_bool(node, "output");
  }
}

void UnstructuredMesh::surface_bins_crossed(
  Position r0, Position r1, const Direction& u, vector<int>& bins) const
{
  fatal_error("Unstructured mesh surface tallies are not implemented.");
}

std::string UnstructuredMesh::bin_label(int bin) const
{
  return fmt::format("Mesh Index ({})", bin);
};

void UnstructuredMesh::to_hdf5(hid_t group) const
{
  hid_t mesh_group = create_group(group, fmt::format("mesh {}", id_));

  write_dataset(mesh_group, "type", "unstructured");
  write_dataset(mesh_group, "filename", filename_);
  write_dataset(mesh_group, "library", this->library());
  // write volume of each element
  vector<double> tet_vols;
  xt::xtensor<double, 2> centroids({static_cast<size_t>(this->n_bins()), 3});
  for (int i = 0; i < this->n_bins(); i++) {
    tet_vols.emplace_back(this->volume(i));
    auto c = this->centroid(i);
    xt::view(centroids, i, xt::all()) = xt::xarray<double>({c.x, c.y, c.z});
  }

  write_dataset(mesh_group, "volumes", tet_vols);
  write_dataset(mesh_group, "centroids", centroids);
  close_group(mesh_group);
}

void StructuredMesh::get_indices(Position r, int* ijk, bool* in_mesh) const
{
  *in_mesh = true;
  for (int i = 0; i < n_dimension_; ++i) {
    ijk[i] = get_index_in_direction(r[i], i);

    if (ijk[i] < 1 || ijk[i] > shape_[i])
      *in_mesh = false;
  }
}

int StructuredMesh::get_bin_from_indices(const int* ijk) const
{
  switch (n_dimension_) {
  case 1:
    return ijk[0] - 1;
  case 2:
    return (ijk[1] - 1) * shape_[0] + ijk[0] - 1;
  case 3:
    return ((ijk[2] - 1) * shape_[1] + (ijk[1] - 1)) * shape_[0] + ijk[0] - 1;
  default:
    throw std::runtime_error {"Invalid number of mesh dimensions"};
  }
}

void StructuredMesh::get_indices_from_bin(int bin, int* ijk) const
{
  if (n_dimension_ == 1) {
    ijk[0] = bin + 1;
  } else if (n_dimension_ == 2) {
    ijk[0] = bin % shape_[0] + 1;
    ijk[1] = bin / shape_[0] + 1;
  } else if (n_dimension_ == 3) {
    ijk[0] = bin % shape_[0] + 1;
    ijk[1] = (bin % (shape_[0] * shape_[1])) / shape_[0] + 1;
    ijk[2] = bin / (shape_[0] * shape_[1]) + 1;
  }
}

int StructuredMesh::get_bin(Position r) const
{
  // Determine indices
  vector<int> ijk(n_dimension_);
  bool in_mesh;
  get_indices(r, ijk.data(), &in_mesh);
  if (!in_mesh)
    return -1;

  // Convert indices to bin
  return get_bin_from_indices(ijk.data());
}

int StructuredMesh::n_bins() const
{
  return xt::prod(shape_)();
}

int StructuredMesh::n_surface_bins() const
{
  return 4 * n_dimension_ * n_bins();
}

xt::xtensor<double, 1> StructuredMesh::count_sites(
  const SourceSite* bank, int64_t length, bool* outside) const
{
  // Determine shape of array for counts
  std::size_t m = this->n_bins();
  vector<std::size_t> shape = {m};

  // Create array of zeros
  xt::xarray<double> cnt {shape, 0.0};
  bool outside_ = false;

  for (int64_t i = 0; i < length; i++) {
    const auto& site = bank[i];

    // determine scoring bin for entropy mesh
    int mesh_bin = get_bin(site.r);

    // if outside mesh, skip particle
    if (mesh_bin < 0) {
      outside_ = true;
      continue;
    }

    // Add to appropriate bin
    cnt(mesh_bin) += site.wgt;
  }

  // Create copy of count data. Since ownership will be acquired by xtensor,
  // std::allocator must be used to avoid Valgrind mismatched free() / delete
  // warnings.
  int total = cnt.size();
  double* cnt_reduced = std::allocator<double> {}.allocate(total);

#ifdef OPENMC_MPI
  // collect values from all processors
  MPI_Reduce(
    cnt.data(), cnt_reduced, total, MPI_DOUBLE, MPI_SUM, 0, mpi::intracomm);

  // Check if there were sites outside the mesh for any processor
  if (outside) {
    MPI_Reduce(&outside_, outside, 1, MPI_C_BOOL, MPI_LOR, 0, mpi::intracomm);
  }
#else
  std::copy(cnt.data(), cnt.data() + total, cnt_reduced);
  if (outside)
    *outside = outside_;
#endif

  // Adapt reduced values in array back into an xarray
  auto arr = xt::adapt(cnt_reduced, total, xt::acquire_ownership(), shape);
  xt::xarray<double> counts = arr;

  return counts;
}

bool StructuredMesh::intersects(Position& r0, Position r1, int* ijk) const
{
  switch (n_dimension_) {
  case 1:
    return intersects_1d(r0, r1, ijk);
  case 2:
    return intersects_2d(r0, r1, ijk);
  case 3:
    return intersects_3d(r0, r1, ijk);
  default:
    throw std::runtime_error {"Invalid number of mesh dimensions."};
  }
}

bool StructuredMesh::intersects_1d(Position& r0, Position r1, int* ijk) const
{
  // Copy coordinates of starting point
  double x0 = r0.x;
  double y0 = r0.y;
  double z0 = r0.z;

  // Copy coordinates of ending point
  double x1 = r1.x;
  double y1 = r1.y;
  double z1 = r1.z;

  // Copy coordinates of mesh lower_left and upper_right
  double xm0 = lower_left_[0];
  double xm1 = upper_right_[0];

  double min_dist = INFTY;

  // Check if line intersects left surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm0 && x1 > xm0) || (x0 > xm0 && x1 < xm0)) {
    double yi = y0 + (xm0 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm0 - x0) * (z1 - z0) / (x1 - x0);
    if (check_intersection_point(xm0, x0, yi, yi, zi, zi, r0, min_dist)) {
      ijk[0] = 1;
    }
  }

  // Check if line intersects right surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm1 && x1 > xm1) || (x0 > xm1 && x1 < xm1)) {
    double yi = y0 + (xm1 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm1 - x0) * (z1 - z0) / (x1 - x0);
    if (check_intersection_point(xm1, x0, yi, yi, zi, zi, r0, min_dist)) {
      ijk[0] = shape_[0];
    }
  }

  return min_dist < INFTY;
}

bool StructuredMesh::intersects_2d(Position& r0, Position r1, int* ijk) const
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

  // Copy coordinates of mesh upper_right
  double xm1 = upper_right_[0];
  double ym1 = upper_right_[1];

  double min_dist = INFTY;

  // Check if line intersects left surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm0 && x1 > xm0) || (x0 > xm0 && x1 < xm0)) {
    double yi = y0 + (xm0 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm0 - x0) * (z1 - z0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1) {
      if (check_intersection_point(xm0, x0, yi, y0, zi, zi, r0, min_dist)) {
        ijk[0] = 1;
        ijk[1] = get_index_in_direction(yi, 1);
      }
    }
  }

  // Check if line intersects back surface -- calculate the intersection point
  // (x,z)
  if ((y0 < ym0 && y1 > ym0) || (y0 > ym0 && y1 < ym0)) {
    double xi = x0 + (ym0 - y0) * (x1 - x0) / (y1 - y0);
    double zi = z0 + (ym0 - y0) * (z1 - z0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1) {
      if (check_intersection_point(xi, x0, ym0, y0, zi, zi, r0, min_dist)) {
        ijk[0] = get_index_in_direction(xi, 0);
        ijk[1] = 1;
      }
    }
  }

  // Check if line intersects right surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm1 && x1 > xm1) || (x0 > xm1 && x1 < xm1)) {
    double yi = y0 + (xm1 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm1 - x0) * (z1 - z0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1) {
      if (check_intersection_point(xm1, x0, yi, y0, zi, zi, r0, min_dist)) {
        ijk[0] = shape_[0];
        ijk[1] = get_index_in_direction(yi, 1);
      }
    }
  }

  // Check if line intersects front surface -- calculate the intersection point
  // (x,z)
  if ((y0 < ym1 && y1 > ym1) || (y0 > ym1 && y1 < ym1)) {
    double xi = x0 + (ym1 - y0) * (x1 - x0) / (y1 - y0);
    double zi = z0 + (ym1 - y0) * (z1 - z0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1) {
      if (check_intersection_point(xi, x0, ym1, y0, zi, zi, r0, min_dist)) {
        ijk[0] = get_index_in_direction(xi, 0);
        ijk[1] = shape_[1];
      }
    }
  }

  return min_dist < INFTY;
}

bool StructuredMesh::intersects_3d(Position& r0, Position r1, int* ijk) const
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

  double min_dist = INFTY;

  // Check if line intersects left surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm0 && x1 > xm0) || (x0 > xm0 && x1 < xm0)) {
    double yi = y0 + (xm0 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm0 - x0) * (z1 - z0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1 && zi >= zm0 && zi < zm1) {
      if (check_intersection_point(xm0, x0, yi, y0, zi, z0, r0, min_dist)) {
        ijk[0] = 1;
        ijk[1] = get_index_in_direction(yi, 1);
        ijk[2] = get_index_in_direction(zi, 2);
      }
    }
  }

  // Check if line intersects back surface -- calculate the intersection point
  // (x,z)
  if ((y0 < ym0 && y1 > ym0) || (y0 > ym0 && y1 < ym0)) {
    double xi = x0 + (ym0 - y0) * (x1 - x0) / (y1 - y0);
    double zi = z0 + (ym0 - y0) * (z1 - z0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1 && zi >= zm0 && zi < zm1) {
      if (check_intersection_point(xi, x0, ym0, y0, zi, z0, r0, min_dist)) {
        ijk[0] = get_index_in_direction(xi, 0);
        ijk[1] = 1;
        ijk[2] = get_index_in_direction(zi, 2);
      }
    }
  }

  // Check if line intersects bottom surface -- calculate the intersection
  // point (x,y)
  if ((z0 < zm0 && z1 > zm0) || (z0 > zm0 && z1 < zm0)) {
    double xi = x0 + (zm0 - z0) * (x1 - x0) / (z1 - z0);
    double yi = y0 + (zm0 - z0) * (y1 - y0) / (z1 - z0);
    if (xi >= xm0 && xi < xm1 && yi >= ym0 && yi < ym1) {
      if (check_intersection_point(xi, x0, yi, y0, zm0, z0, r0, min_dist)) {
        ijk[0] = get_index_in_direction(xi, 0);
        ijk[1] = get_index_in_direction(yi, 1);
        ijk[2] = 1;
      }
    }
  }

  // Check if line intersects right surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm1 && x1 > xm1) || (x0 > xm1 && x1 < xm1)) {
    double yi = y0 + (xm1 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm1 - x0) * (z1 - z0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1 && zi >= zm0 && zi < zm1) {
      if (check_intersection_point(xm1, x0, yi, y0, zi, z0, r0, min_dist)) {
        ijk[0] = shape_[0];
        ijk[1] = get_index_in_direction(yi, 1);
        ijk[2] = get_index_in_direction(zi, 2);
      }
    }
  }

  // Check if line intersects front surface -- calculate the intersection point
  // (x,z)
  if ((y0 < ym1 && y1 > ym1) || (y0 > ym1 && y1 < ym1)) {
    double xi = x0 + (ym1 - y0) * (x1 - x0) / (y1 - y0);
    double zi = z0 + (ym1 - y0) * (z1 - z0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1 && zi >= zm0 && zi < zm1) {
      if (check_intersection_point(xi, x0, ym1, y0, zi, z0, r0, min_dist)) {
        ijk[0] = get_index_in_direction(xi, 0);
        ijk[1] = shape_[1];
        ijk[2] = get_index_in_direction(zi, 2);
      }
    }
  }

  // Check if line intersects top surface -- calculate the intersection point
  // (x,y)
  if ((z0 < zm1 && z1 > zm1) || (z0 > zm1 && z1 < zm1)) {
    double xi = x0 + (zm1 - z0) * (x1 - x0) / (z1 - z0);
    double yi = y0 + (zm1 - z0) * (y1 - y0) / (z1 - z0);
    if (xi >= xm0 && xi < xm1 && yi >= ym0 && yi < ym1) {
      if (check_intersection_point(xi, x0, yi, y0, zm1, z0, r0, min_dist)) {
        ijk[0] = get_index_in_direction(xi, 0);
        ijk[1] = get_index_in_direction(yi, 1);
        ijk[2] = shape_[2];
      }
    }
  }

  return min_dist < INFTY;
}

void StructuredMesh::bins_crossed(Position r0, Position r1, const Direction& u,
  vector<int>& bins, vector<double>& lengths) const
{
  // ========================================================================
  // Determine where the track intersects the mesh and if it intersects at all.

  // Compute the length of the entire track.
  double total_distance = (r1 - r0).norm();
  if (total_distance == 0.0)
    return;

  // While determining if this track intersects the mesh, offset the starting
  // and ending coords by a bit.  This avoid finite-precision errors that can
  // occur when the mesh surfaces coincide with lattice or geometric surfaces.
  Position last_r = r0 + TINY_BIT * u;
  Position r = r1 - TINY_BIT * u;

  // Determine the mesh indices for the starting and ending coords. Here, we
  // use arrays for ijk0 and ijk1 instead of vector because we obtain a
  // small performance improvement by forcing this data to live on the stack,
  // rather than on the heap. We know the maximum length is 3, and by
  // ensuring that all loops are only indexed up to n_dimension, we will not
  // access any non-initialized values. The same concept is used throughout.
  int n = n_dimension_;
  int ijk0[3], ijk1[3];
  bool start_in_mesh;
  get_indices(last_r, ijk0, &start_in_mesh);
  bool end_in_mesh;
  get_indices(r, ijk1, &end_in_mesh);

  // Reset coordinates and check for a mesh intersection if necessary.
  if (start_in_mesh) {
    // The initial coords lie in the mesh, use those coords for tallying.
    last_r = r0;
  } else {
    // The initial coords do not lie in the mesh.  Check to see if the particle
    // eventually intersects the mesh and compute the relevant coords and
    // indices.
    if (!intersects(last_r, r, ijk0))
      return;
  }
  r = r1;

  // The TINY_BIT offsets above mean that the preceding logic cannot always find
  // the correct ijk0 and ijk1 indices. For tracks shorter than 2*TINY_BIT, just
  // assume the track lies in only one mesh bin. These tracks are very short so
  // any error caused by this assumption will be small. It is important that
  // ijk0 values are used rather than ijk1 because the previous logic guarantees
  // ijk0 is a valid mesh bin.
  if (total_distance < 2 * TINY_BIT) {
    for (int i = 0; i < n; ++i)
      ijk1[i] = ijk0[i];
  }

  // ========================================================================
  // Find which mesh cells are traversed and the length of each traversal.

  while (true) {
    if (std::equal(ijk0, ijk0 + n, ijk1)) {
      // The track ends in this cell.  Use the particle end location rather
      // than the mesh surface and stop iterating.
      double distance = (last_r - r).norm();
      bins.push_back(get_bin_from_indices(ijk0));
      lengths.push_back(distance / total_distance);
      break;
    }

    // The track exits this cell.  Determine the distance to each mesh surface.
    double d[3];
    for (int k = 0; k < n; ++k) {
      if (std::fabs(u[k]) < FP_PRECISION) {
        d[k] = INFTY;
      } else if (u[k] > 0) {
        double xyz_cross = positive_grid_boundary(ijk0, k);
        d[k] = (xyz_cross - last_r[k]) / u[k];
      } else {
        double xyz_cross = negative_grid_boundary(ijk0, k);
        d[k] = (xyz_cross - last_r[k]) / u[k];
      }
    }

    // Pick the closest mesh surface and append this traversal to the output.
    auto j = std::min_element(d, d + n) - d;
    double distance = d[j];
    bins.push_back(get_bin_from_indices(ijk0));
    lengths.push_back(distance / total_distance);

    // Translate to the oncoming mesh surface.
    last_r += distance * u;

    // Increment the indices into the next mesh cell.
    if (u[j] > 0.0) {
      ++ijk0[j];
    } else {
      --ijk0[j];
    }

    // If the next indices are invalid, then the track has left the mesh and
    // we are done.
    bool in_mesh = true;
    for (int i = 0; i < n; ++i) {
      if (ijk0[i] < 1 || ijk0[i] > shape_[i]) {
        in_mesh = false;
        break;
      }
    }
    if (!in_mesh)
      break;
  }
}

//==============================================================================
// RegularMesh implementation
//==============================================================================

RegularMesh::RegularMesh(pugi::xml_node node) : StructuredMesh {node}
{
  // Determine number of dimensions for mesh
  if (!check_for_node(node, "dimension")) {
    fatal_error("Must specify <dimension> on a regular mesh.");
  }

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

  // Check for lower-left coordinates
  if (check_for_node(node, "lower_left")) {
    // Read mesh lower-left corner location
    lower_left_ = get_node_xarray<double>(node, "lower_left");
  } else {
    fatal_error("Must specify <lower_left> on a mesh.");
  }

  // Make sure lower_left and dimension match
  if (shape_.size() != lower_left_.size()) {
    fatal_error("Number of entries on <lower_left> must be the same "
                "as the number of entries on <dimension>.");
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

    // Set width
    width_ = xt::eval((upper_right_ - lower_left_) / shape_);
  } else {
    fatal_error("Must specify either <upper_right> and <width> on a mesh.");
  }

  // Set volume fraction
  volume_frac_ = 1.0 / xt::prod(shape_)();
}

int RegularMesh::get_index_in_direction(double r, int i) const
{
  return std::ceil((r - lower_left_[i]) / width_[i]);
}

double RegularMesh::positive_grid_boundary(int* ijk, int i) const
{
  return lower_left_[i] + ijk[i] * width_[i];
}

double RegularMesh::negative_grid_boundary(int* ijk, int i) const
{
  return lower_left_[i] + (ijk[i] - 1) * width_[i];
}

void RegularMesh::surface_bins_crossed(
  Position r0, Position r1, const Direction& u, vector<int>& bins) const
{
  // Determine indices for starting and ending location.
  int n = n_dimension_;
  vector<int> ijk0(n), ijk1(n);
  bool start_in_mesh;
  get_indices(r0, ijk0.data(), &start_in_mesh);
  bool end_in_mesh;
  get_indices(r1, ijk1.data(), &end_in_mesh);

  // Check if the track intersects any part of the mesh.
  if (!start_in_mesh) {
    Position r0_copy = r0;
    vector<int> ijk0_copy(ijk0);
    if (!intersects(r0_copy, r1, ijk0_copy.data()))
      return;
  }

  // ========================================================================
  // Find which mesh surfaces are crossed.

  // Calculate number of surface crossings
  int n_cross = 0;
  for (int i = 0; i < n; ++i)
    n_cross += std::abs(ijk1[i] - ijk0[i]);
  if (n_cross == 0)
    return;

  // Bounding coordinates
  Position xyz_cross;
  for (int i = 0; i < n; ++i) {
    if (u[i] > 0.0) {
      xyz_cross[i] = positive_grid_boundary(ijk0.data(), i);
    } else {
      xyz_cross[i] = negative_grid_boundary(ijk0.data(), i);
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
        d[i] = INFTY;
      } else {
        d[i] = (xyz_cross[i] - r0[i]) / u[i];
      }
      distance = std::min(distance, d[i]);
    }

    // Loop over the dimensions
    for (int i = 0; i < n; ++i) {
      // Check whether distance is the shortest distance
      if (distance == d[i]) {

        // Check whether the current indices are within the mesh bounds
        bool in_mesh = true;
        for (int j = 0; j < n; ++j) {
          if (ijk0[j] < 1 || ijk0[j] > shape_[j]) {
            in_mesh = false;
            break;
          }
        }

        // Check whether particle is moving in positive i direction
        if (u[i] > 0) {

          // Outward current on i max surface
          if (in_mesh) {
            int i_surf = 4 * i + 3;
            int i_mesh = get_bin_from_indices(ijk0.data());
            int i_bin = 4 * n * i_mesh + i_surf - 1;

            bins.push_back(i_bin);
          }

          // Advance position
          ++ijk0[i];
          xyz_cross[i] += width_[i];
          in_mesh = true;
          for (int j = 0; j < n; ++j) {
            if (ijk0[j] < 1 || ijk0[j] > shape_[j]) {
              in_mesh = false;
              break;
            }
          }

          // If the particle crossed the surface, tally the inward current on
          // i min surface
          if (in_mesh) {
            int i_surf = 4 * i + 2;
            int i_mesh = get_bin_from_indices(ijk0.data());
            int i_bin = 4 * n * i_mesh + i_surf - 1;

            bins.push_back(i_bin);
          }

        } else {
          // The particle is moving in the negative i direction

          // Outward current on i min surface
          if (in_mesh) {
            int i_surf = 4 * i + 1;
            int i_mesh = get_bin_from_indices(ijk0.data());
            int i_bin = 4 * n * i_mesh + i_surf - 1;

            bins.push_back(i_bin);
          }

          // Advance position
          --ijk0[i];
          xyz_cross[i] -= width_[i];
          in_mesh = true;
          for (int j = 0; j < n; ++j) {
            if (ijk0[j] < 1 || ijk0[j] > shape_[j]) {
              in_mesh = false;
              break;
            }
          }

          // If the particle crossed the surface, tally the inward current on
          // i max surface
          if (in_mesh) {
            int i_surf = 4 * i + 4;
            int i_mesh = get_bin_from_indices(ijk0.data());
            int i_bin = 4 * n * i_mesh + i_surf - 1;

            bins.push_back(i_bin);
          }
        }
      }
    }

    // Calculate new coordinates
    r0 += distance * u;
  }
}

std::pair<vector<double>, vector<double>> RegularMesh::plot(
  Position plot_ll, Position plot_ur) const
{
  // Figure out which axes lie in the plane of the plot.
  array<int, 2> axes {-1, -1};
  if (plot_ur.z == plot_ll.z) {
    axes[0] = 0;
    if (n_dimension_ > 1)
      axes[1] = 1;
  } else if (plot_ur.y == plot_ll.y) {
    axes[0] = 0;
    if (n_dimension_ > 2)
      axes[1] = 2;
  } else if (plot_ur.x == plot_ll.x) {
    if (n_dimension_ > 1)
      axes[0] = 1;
    if (n_dimension_ > 2)
      axes[1] = 2;
  } else {
    fatal_error("Can only plot mesh lines on an axis-aligned plot");
  }

  // Get the coordinates of the mesh lines along both of the axes.
  array<vector<double>, 2> axis_lines;
  for (int i_ax = 0; i_ax < 2; ++i_ax) {
    int axis = axes[i_ax];
    if (axis == -1)
      continue;
    auto& lines {axis_lines[i_ax]};

    double coord = lower_left_[axis];
    for (int i = 0; i < shape_[axis] + 1; ++i) {
      if (coord >= plot_ll[axis] && coord <= plot_ur[axis])
        lines.push_back(coord);
      coord += width_[axis];
    }
  }

  return {axis_lines[0], axis_lines[1]};
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

xt::xtensor<double, 1> RegularMesh::count_sites(
  const SourceSite* bank, int64_t length, bool* outside) const
{
  // Determine shape of array for counts
  std::size_t m = this->n_bins();
  vector<std::size_t> shape = {m};

  // Create array of zeros
  xt::xarray<double> cnt {shape, 0.0};
  bool outside_ = false;

  for (int64_t i = 0; i < length; i++) {
    const auto& site = bank[i];

    // determine scoring bin for entropy mesh
    int mesh_bin = get_bin(site.r);

    // if outside mesh, skip particle
    if (mesh_bin < 0) {
      outside_ = true;
      continue;
    }

    // Add to appropriate bin
    cnt(mesh_bin) += site.wgt;
  }

  // Create copy of count data. Since ownership will be acquired by xtensor,
  // std::allocator must be used to avoid Valgrind mismatched free() / delete
  // warnings.
  int total = cnt.size();
  double* cnt_reduced = std::allocator<double> {}.allocate(total);

#ifdef OPENMC_MPI
  // collect values from all processors
  MPI_Reduce(
    cnt.data(), cnt_reduced, total, MPI_DOUBLE, MPI_SUM, 0, mpi::intracomm);

  // Check if there were sites outside the mesh for any processor
  if (outside) {
    MPI_Reduce(&outside_, outside, 1, MPI_C_BOOL, MPI_LOR, 0, mpi::intracomm);
  }
#else
  std::copy(cnt.data(), cnt.data() + total, cnt_reduced);
  if (outside)
    *outside = outside_;
#endif

  // Adapt reduced values in array back into an xarray
  auto arr = xt::adapt(cnt_reduced, total, xt::acquire_ownership(), shape);
  xt::xarray<double> counts = arr;

  return counts;
}

//==============================================================================
// RectilinearMesh implementation
//==============================================================================

RectilinearMesh::RectilinearMesh(pugi::xml_node node) : StructuredMesh {node}
{
  n_dimension_ = 3;

  grid_.resize(n_dimension_);
  grid_[0] = get_node_array<double>(node, "x_grid");
  grid_[1] = get_node_array<double>(node, "y_grid");
  grid_[2] = get_node_array<double>(node, "z_grid");

  if (int err = set_grid()) {
    fatal_error(openmc_err_msg);
  }
}

double RectilinearMesh::positive_grid_boundary(int* ijk, int i) const
{
  return grid_[i][ijk[i]];
}

double RectilinearMesh::negative_grid_boundary(int* ijk, int i) const
{
  return grid_[i][ijk[i] - 1];
}

int RectilinearMesh::set_grid()
{
  shape_ = {static_cast<int>(grid_[0].size()) - 1,
    static_cast<int>(grid_[1].size()) - 1,
    static_cast<int>(grid_[2].size()) - 1};

  for (const auto& g : grid_) {
    if (g.size() < 2) {
      set_errmsg("x-, y-, and z- grids for rectilinear meshes "
                 "must each have at least 2 points");
      return OPENMC_E_INVALID_ARGUMENT;
    }
    for (int i = 1; i < g.size(); ++i) {
      if (g[i] <= g[i - 1]) {
        set_errmsg("Values in for x-, y-, and z- grids for "
                   "rectilinear meshes must be sorted and unique.");
        return OPENMC_E_INVALID_ARGUMENT;
      }
    }
  }

  lower_left_ = {grid_[0].front(), grid_[1].front(), grid_[2].front()};
  upper_right_ = {grid_[0].back(), grid_[1].back(), grid_[2].back()};

  return 0;
}

void RectilinearMesh::surface_bins_crossed(
  Position r0, Position r1, const Direction& u, vector<int>& bins) const
{
  // Determine indices for starting and ending location.
  int ijk0[3], ijk1[3];
  bool start_in_mesh;
  get_indices(r0, ijk0, &start_in_mesh);
  bool end_in_mesh;
  get_indices(r1, ijk1, &end_in_mesh);

  // If the starting coordinates do not lie in the mesh, compute the coords and
  // mesh indices of the first intersection, and add the bin for this first
  // intersection.  Return if the particle does not intersect the mesh at all.
  if (!start_in_mesh) {
    // Compute the incoming intersection coordinates and indices.
    if (!intersects(r0, r1, ijk0))
      return;

    // Determine which surface the particle entered.
    double min_dist = INFTY;
    int i_surf;
    for (int i = 0; i < 3; ++i) {
      if (u[i] > 0.0 && ijk0[i] == 1) {
        double d = std::abs(r0[i] - grid_[i][0]);
        if (d < min_dist) {
          min_dist = d;
          i_surf = 4 * i + 2;
        }
      } else if (u[i] < 0.0 && ijk0[i] == shape_[i]) {
        double d = std::abs(r0[i] - grid_[i][shape_[i]]);
        if (d < min_dist) {
          min_dist = d;
          i_surf = 4 * i + 4;
        }
      } // u[i] == 0 intentionally skipped
    }

    // Add the incoming current bin.
    int i_mesh = get_bin_from_indices(ijk0);
    int i_bin = 4 * 3 * i_mesh + i_surf - 1;
    bins.push_back(i_bin);
  }

  // If the ending coordinates do not lie in the mesh, compute the coords and
  // mesh indices of the last intersection, and add the bin for this last
  // intersection.
  if (!end_in_mesh) {
    // Compute the outgoing intersection coordinates and indices.
    intersects(r1, r0, ijk1);

    // Determine which surface the particle exited.
    double min_dist = INFTY;
    int i_surf;
    for (int i = 0; i < 3; ++i) {
      if (u[i] > 0.0 && ijk1[i] == shape_[i]) {
        double d = std::abs(r1[i] - grid_[i][shape_[i]]);
        if (d < min_dist) {
          min_dist = d;
          i_surf = 4 * i + 3;
        }
      } else if (u[i] < 0.0 && ijk1[i] == 1) {
        double d = std::abs(r1[i] - grid_[i][0]);
        if (d < min_dist) {
          min_dist = d;
          i_surf = 4 * i + 1;
        }
      } // u[i] == 0 intentionally skipped
    }

    // Add the outgoing current bin.
    int i_mesh = get_bin_from_indices(ijk1);
    int i_bin = 4 * 3 * i_mesh + i_surf - 1;
    bins.push_back(i_bin);
  }

  // ========================================================================
  // Find which mesh surfaces are crossed.

  // Calculate number of surface crossings
  int n_cross = 0;
  for (int i = 0; i < 3; ++i)
    n_cross += std::abs(ijk1[i] - ijk0[i]);
  if (n_cross == 0)
    return;

  // Bounding coordinates
  Position xyz_cross;
  for (int i = 0; i < 3; ++i) {
    if (u[i] > 0.0) {
      xyz_cross[i] = positive_grid_boundary(ijk0, i);
    } else {
      xyz_cross[i] = negative_grid_boundary(ijk0, i);
    }
  }

  for (int j = 0; j < n_cross; ++j) {
    // Set the distances to infinity
    Position d {INFTY, INFTY, INFTY};

    // Determine closest bounding surface. We need to treat
    // special case where the cosine of the angle is zero since this would
    // result in a divide-by-zero.
    double distance = INFTY;
    for (int i = 0; i < 3; ++i) {
      if (u[i] == 0) {
        d[i] = INFTY;
      } else {
        d[i] = (xyz_cross[i] - r0[i]) / u[i];
      }
      distance = std::min(distance, d[i]);
    }

    // Loop over the dimensions
    for (int i = 0; i < 3; ++i) {
      // Check whether distance is the shortest distance
      if (distance == d[i]) {

        // Check whether particle is moving in positive i direction
        if (u[i] > 0) {

          // Outward current on i max surface
          int i_surf = 4 * i + 3;
          int i_mesh = get_bin_from_indices(ijk0);
          int i_bin = 4 * 3 * i_mesh + i_surf - 1;
          bins.push_back(i_bin);

          // Advance position
          ++ijk0[i];
          xyz_cross[i] = grid_[i][ijk0[i]];

          // Inward current on i min surface
          i_surf = 4 * i + 2;
          i_mesh = get_bin_from_indices(ijk0);
          i_bin = 4 * 3 * i_mesh + i_surf - 1;
          bins.push_back(i_bin);

        } else {
          // The particle is moving in the negative i direction

          // Outward current on i min surface
          int i_surf = 4 * i + 1;
          int i_mesh = get_bin_from_indices(ijk0);
          int i_bin = 4 * 3 * i_mesh + i_surf - 1;
          bins.push_back(i_bin);

          // Advance position
          --ijk0[i];
          xyz_cross[i] = grid_[i][ijk0[i] - 1];

          // Inward current on i min surface
          i_surf = 4 * i + 4;
          i_mesh = get_bin_from_indices(ijk0);
          i_bin = 4 * 3 * i_mesh + i_surf - 1;
          bins.push_back(i_bin);
        }
      }
    }

    // Calculate new coordinates
    r0 += distance * u;
  }
}

int RectilinearMesh::get_index_in_direction(double r, int i) const
{
  return lower_bound_index(grid_[i].begin(), grid_[i].end(), r) + 1;
}

std::pair<vector<double>, vector<double>> RectilinearMesh::plot(
  Position plot_ll, Position plot_ur) const
{
  // Figure out which axes lie in the plane of the plot.
  array<int, 2> axes {-1, -1};
  if (plot_ur.z == plot_ll.z) {
    axes = {0, 1};
  } else if (plot_ur.y == plot_ll.y) {
    axes = {0, 2};
  } else if (plot_ur.x == plot_ll.x) {
    axes = {1, 2};
  } else {
    fatal_error("Can only plot mesh lines on an axis-aligned plot");
  }

  // Get the coordinates of the mesh lines along both of the axes.
  array<vector<double>, 2> axis_lines;
  for (int i_ax = 0; i_ax < 2; ++i_ax) {
    int axis = axes[i_ax];
    vector<double>& lines {axis_lines[i_ax]};

    for (auto coord : grid_[axis]) {
      if (coord >= plot_ll[axis] && coord <= plot_ur[axis])
        lines.push_back(coord);
    }
  }

  return {axis_lines[0], axis_lines[1]};
}

void RectilinearMesh::to_hdf5(hid_t group) const
{
  hid_t mesh_group = create_group(group, "mesh " + std::to_string(id_));

  write_dataset(mesh_group, "type", "rectilinear");
  write_dataset(mesh_group, "x_grid", grid_[0]);
  write_dataset(mesh_group, "y_grid", grid_[1]);
  write_dataset(mesh_group, "z_grid", grid_[2]);

  close_group(mesh_group);
}

//==============================================================================
// Helper functions for the C API
//==============================================================================

int check_mesh(int32_t index)
{
  if (index < 0 || index >= model::meshes.size()) {
    set_errmsg("Index in meshes array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  return 0;
}

template<class T>
int check_mesh_type(int32_t index)
{
  if (int err = check_mesh(index))
    return err;

  T* mesh = dynamic_cast<T*>(model::meshes[index].get());
  if (!mesh) {
    set_errmsg("This function is not valid for input mesh.");
    return OPENMC_E_INVALID_TYPE;
  }
  return 0;
}

//==============================================================================
// C API functions
//==============================================================================

// Return the type of mesh as a C string
extern "C" int openmc_mesh_get_type(int32_t index, char* type)
{
  if (int err = check_mesh(index))
    return err;

  RegularMesh* mesh = dynamic_cast<RegularMesh*>(model::meshes[index].get());
  if (mesh) {
    std::strcpy(type, "regular");
  } else {
    RectilinearMesh* mesh =
      dynamic_cast<RectilinearMesh*>(model::meshes[index].get());
    if (mesh) {
      std::strcpy(type, "rectilinear");
    }
  }
  return 0;
}

//! Extend the meshes array by n elements
extern "C" int openmc_extend_meshes(
  int32_t n, const char* type, int32_t* index_start, int32_t* index_end)
{
  if (index_start)
    *index_start = model::meshes.size();
  std::string mesh_type;

  for (int i = 0; i < n; ++i) {
    if (std::strcmp(type, "regular") == 0) {
      model::meshes.push_back(make_unique<RegularMesh>());
    } else if (std::strcmp(type, "rectilinear") == 0) {
      model::meshes.push_back(make_unique<RectilinearMesh>());
    } else {
      throw std::runtime_error {"Unknown mesh type: " + std::string(type)};
    }
  }
  if (index_end)
    *index_end = model::meshes.size() - 1;

  return 0;
}

//! Adds a new unstructured mesh to OpenMC
extern "C" int openmc_add_unstructured_mesh(
  const char filename[], const char library[], int* id)
{
  std::string lib_name(library);
  std::string mesh_file(filename);
  bool valid_lib = false;

#ifdef DAGMC
  if (lib_name == "moab") {
    model::meshes.push_back(std::move(make_unique<MOABMesh>(mesh_file)));
    valid_lib = true;
  }
#endif

#ifdef LIBMESH
  if (lib_name == "libmesh") {
    model::meshes.push_back(std::move(make_unique<LibMesh>(mesh_file)));
    valid_lib = true;
  }
#endif

  if (!valid_lib) {
    set_errmsg(fmt::format("Mesh library {} is not supported "
                           "by this build of OpenMC",
      lib_name));
    return OPENMC_E_INVALID_ARGUMENT;
  }

  // auto-assign new ID
  model::meshes.back()->set_id(-1);
  *id = model::meshes.back()->id_;

  return 0;
}

//! Return the index in the meshes array of a mesh with a given ID
extern "C" int openmc_get_mesh_index(int32_t id, int32_t* index)
{
  auto pair = model::mesh_map.find(id);
  if (pair == model::mesh_map.end()) {
    set_errmsg("No mesh exists with ID=" + std::to_string(id) + ".");
    return OPENMC_E_INVALID_ID;
  }
  *index = pair->second;
  return 0;
}

//! Return the ID of a mesh
extern "C" int openmc_mesh_get_id(int32_t index, int32_t* id)
{
  if (int err = check_mesh(index))
    return err;
  *id = model::meshes[index]->id_;
  return 0;
}

//! Set the ID of a mesh
extern "C" int openmc_mesh_set_id(int32_t index, int32_t id)
{
  if (int err = check_mesh(index))
    return err;
  model::meshes[index]->id_ = id;
  model::mesh_map[id] = index;
  return 0;
}

//! Get the dimension of a regular mesh
extern "C" int openmc_regular_mesh_get_dimension(
  int32_t index, int** dims, int* n)
{
  if (int err = check_mesh_type<RegularMesh>(index))
    return err;
  RegularMesh* mesh = dynamic_cast<RegularMesh*>(model::meshes[index].get());
  *dims = mesh->shape_.data();
  *n = mesh->n_dimension_;
  return 0;
}

//! Set the dimension of a regular mesh
extern "C" int openmc_regular_mesh_set_dimension(
  int32_t index, int n, const int* dims)
{
  if (int err = check_mesh_type<RegularMesh>(index))
    return err;
  RegularMesh* mesh = dynamic_cast<RegularMesh*>(model::meshes[index].get());

  // Copy dimension
  vector<std::size_t> shape = {static_cast<std::size_t>(n)};
  mesh->shape_ = xt::adapt(dims, n, xt::no_ownership(), shape);
  mesh->n_dimension_ = mesh->shape_.size();
  return 0;
}

//! Get the regular mesh parameters
extern "C" int openmc_regular_mesh_get_params(
  int32_t index, double** ll, double** ur, double** width, int* n)
{
  if (int err = check_mesh_type<RegularMesh>(index))
    return err;
  RegularMesh* m = dynamic_cast<RegularMesh*>(model::meshes[index].get());

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

//! Set the regular mesh parameters
extern "C" int openmc_regular_mesh_set_params(
  int32_t index, int n, const double* ll, const double* ur, const double* width)
{
  if (int err = check_mesh_type<RegularMesh>(index))
    return err;
  RegularMesh* m = dynamic_cast<RegularMesh*>(model::meshes[index].get());

  vector<std::size_t> shape = {static_cast<std::size_t>(n)};
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

//! Get the rectilinear mesh grid
extern "C" int openmc_rectilinear_mesh_get_grid(int32_t index, double** grid_x,
  int* nx, double** grid_y, int* ny, double** grid_z, int* nz)
{
  if (int err = check_mesh_type<RectilinearMesh>(index))
    return err;
  RectilinearMesh* m =
    dynamic_cast<RectilinearMesh*>(model::meshes[index].get());

  if (m->lower_left_.dimension() == 0) {
    set_errmsg("Mesh parameters have not been set.");
    return OPENMC_E_ALLOCATE;
  }

  *grid_x = m->grid_[0].data();
  *nx = m->grid_[0].size();
  *grid_y = m->grid_[1].data();
  *ny = m->grid_[1].size();
  *grid_z = m->grid_[2].data();
  *nz = m->grid_[2].size();

  return 0;
}

//! Set the rectilienar mesh parameters
extern "C" int openmc_rectilinear_mesh_set_grid(int32_t index,
  const double* grid_x, const int nx, const double* grid_y, const int ny,
  const double* grid_z, const int nz)
{
  if (int err = check_mesh_type<RectilinearMesh>(index))
    return err;
  RectilinearMesh* m =
    dynamic_cast<RectilinearMesh*>(model::meshes[index].get());

  m->n_dimension_ = 3;
  m->grid_.resize(m->n_dimension_);

  for (int i = 0; i < nx; i++) {
    m->grid_[0].push_back(grid_x[i]);
  }
  for (int i = 0; i < ny; i++) {
    m->grid_[1].push_back(grid_y[i]);
  }
  for (int i = 0; i < nz; i++) {
    m->grid_[2].push_back(grid_z[i]);
  }

  int err = m->set_grid();
  return err;
}

#ifdef DAGMC

MOABMesh::MOABMesh(pugi::xml_node node) : UnstructuredMesh(node)
{
  initialize();
}

MOABMesh::MOABMesh(const std::string& filename)
{
  filename_ = filename;
  initialize();
}

MOABMesh::MOABMesh(std::shared_ptr<moab::Interface> external_mbi)
{
  mbi_ = external_mbi;
  filename_ = "unknown (external file)";
  this->initialize();
}

void MOABMesh::initialize()
{

  // Create the MOAB interface and load data from file
  this->create_interface();

  // Initialise MOAB error code
  moab::ErrorCode rval = moab::MB_SUCCESS;

  // Set the dimension
  n_dimension_ = 3;

  // set member range of tetrahedral entities
  rval = mbi_->get_entities_by_dimension(0, n_dimension_, ehs_);
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to get all tetrahedral elements");
  }

  if (!ehs_.all_of_type(moab::MBTET)) {
    warning("Non-tetrahedral elements found in unstructured "
            "mesh file: " +
            filename_);
  }

  // make an entity set for all tetrahedra
  // this is used for convenience later in output
  rval = mbi_->create_meshset(moab::MESHSET_SET, tetset_);
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to create an entity set for the tetrahedral elements");
  }

  rval = mbi_->add_entities(tetset_, ehs_);
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to add tetrahedra to an entity set.");
  }

  // build acceleration data structures
  compute_barycentric_data(ehs_);
  build_kdtree(ehs_);
}

void MOABMesh::create_interface()
{
  // Do not create a MOAB instance if one is already in memory
  if (mbi_)
    return;

  // create MOAB instance
  mbi_ = std::make_shared<moab::Core>();

  // load unstructured mesh file
  moab::ErrorCode rval = mbi_->load_file(filename_.c_str());
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to load the unstructured mesh file: " + filename_);
  }
}

void MOABMesh::build_kdtree(const moab::Range& all_tets)
{
  moab::Range all_tris;
  int adj_dim = 2;
  moab::ErrorCode rval = mbi_->get_adjacencies(
    all_tets, adj_dim, true, all_tris, moab::Interface::UNION);
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to get adjacent triangles for tets");
  }

  if (!all_tris.all_of_type(moab::MBTRI)) {
    warning("Non-triangle elements found in tet adjacencies in "
            "unstructured mesh file: " +
            filename_);
  }

  // combine into one range
  moab::Range all_tets_and_tris;
  all_tets_and_tris.merge(all_tets);
  all_tets_and_tris.merge(all_tris);

  // create a kd-tree instance
  kdtree_ = make_unique<moab::AdaptiveKDTree>(mbi_.get());

  // build the tree
  rval = kdtree_->build_tree(all_tets_and_tris, &kdtree_root_);
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to construct KDTree for the "
                "unstructured mesh file: " +
                filename_);
  }
}

void MOABMesh::intersect_track(const moab::CartVect& start,
  const moab::CartVect& dir, double track_len, vector<double>& hits) const
{
  hits.clear();

  moab::ErrorCode rval;
  vector<moab::EntityHandle> tris;
  // get all intersections with triangles in the tet mesh
  // (distances are relative to the start point, not the previous intersection)
  rval = kdtree_->ray_intersect_triangles(kdtree_root_, FP_COINCIDENT,
    dir.array(), start.array(), tris, hits, 0, track_len);
  if (rval != moab::MB_SUCCESS) {
    fatal_error(
      "Failed to compute intersections on unstructured mesh: " + filename_);
  }

  // remove duplicate intersection distances
  std::unique(hits.begin(), hits.end());

  // sorts by first component of std::pair by default
  std::sort(hits.begin(), hits.end());
}

void MOABMesh::bins_crossed(Position r0, Position r1, const Direction& u,
  vector<int>& bins, vector<double>& lengths) const
{
  moab::CartVect start(r0.x, r0.y, r0.z);
  moab::CartVect end(r1.x, r1.y, r1.z);
  moab::CartVect dir(u.x, u.y, u.z);
  dir.normalize();

  double track_len = (end - start).length();
  if (track_len == 0.0)
    return;

  start -= TINY_BIT * dir;
  end += TINY_BIT * dir;

  vector<double> hits;
  intersect_track(start, dir, track_len, hits);

  bins.clear();
  lengths.clear();

  // if there are no intersections the track may lie entirely
  // within a single tet. If this is the case, apply entire
  // score to that tet and return.
  if (hits.size() == 0) {
    Position midpoint = r0 + u * (track_len * 0.5);
    int bin = this->get_bin(midpoint);
    if (bin != -1) {
      bins.push_back(bin);
      lengths.push_back(1.0);
    }
    return;
  }

  // for each segment in the set of tracks, try to look up a tet
  // at the midpoint of the segment
  Position current = r0;
  double last_dist = 0.0;
  for (const auto& hit : hits) {
    // get the segment length
    double segment_length = hit - last_dist;
    last_dist = hit;
    // find the midpoint of this segment
    Position midpoint = current + u * (segment_length * 0.5);
    // try to find a tet for this position
    int bin = this->get_bin(midpoint);

    // determine the start point for this segment
    current = r0 + u * hit;

    if (bin == -1) {
      continue;
    }

    bins.push_back(bin);
    lengths.push_back(segment_length / track_len);
  }

  // tally remaining portion of track after last hit if
  // the last segment of the track is in the mesh but doesn't
  // reach the other side of the tet
  if (hits.back() < track_len) {
    Position segment_start = r0 + u * hits.back();
    double segment_length = track_len - hits.back();
    Position midpoint = segment_start + u * (segment_length * 0.5);
    int bin = this->get_bin(midpoint);
    if (bin != -1) {
      bins.push_back(bin);
      lengths.push_back(segment_length / track_len);
    }
  }
};

moab::EntityHandle MOABMesh::get_tet(const Position& r) const
{
  moab::CartVect pos(r.x, r.y, r.z);
  // find the leaf of the kd-tree for this position
  moab::AdaptiveKDTreeIter kdtree_iter;
  moab::ErrorCode rval = kdtree_->point_search(pos.array(), kdtree_iter);
  if (rval != moab::MB_SUCCESS) {
    return 0;
  }

  // retrieve the tet elements of this leaf
  moab::EntityHandle leaf = kdtree_iter.handle();
  moab::Range tets;
  rval = mbi_->get_entities_by_dimension(leaf, 3, tets, false);
  if (rval != moab::MB_SUCCESS) {
    warning("MOAB error finding tets.");
  }

  // loop over the tets in this leaf, returning the containing tet if found
  for (const auto& tet : tets) {
    if (point_in_tet(pos, tet)) {
      return tet;
    }
  }

  // if no tet is found, return an invalid handle
  return 0;
}

double MOABMesh::volume(int bin) const
{
  return tet_volume(get_ent_handle_from_bin(bin));
}

std::string MOABMesh::library() const
{
  return "moab";
}

double MOABMesh::tet_volume(moab::EntityHandle tet) const
{
  vector<moab::EntityHandle> conn;
  moab::ErrorCode rval = mbi_->get_connectivity(&tet, 1, conn);
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to get tet connectivity");
  }

  moab::CartVect p[4];
  rval = mbi_->get_coords(conn.data(), conn.size(), p[0].array());
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to get tet coords");
  }

  return 1.0 / 6.0 * (((p[1] - p[0]) * (p[2] - p[0])) % (p[3] - p[0]));
}

int MOABMesh::get_bin(Position r) const
{
  moab::EntityHandle tet = get_tet(r);
  if (tet == 0) {
    return -1;
  } else {
    return get_bin_from_ent_handle(tet);
  }
}

void MOABMesh::compute_barycentric_data(const moab::Range& tets)
{
  moab::ErrorCode rval;

  baryc_data_.clear();
  baryc_data_.resize(tets.size());

  // compute the barycentric data for each tet element
  // and store it as a 3x3 matrix
  for (auto& tet : tets) {
    vector<moab::EntityHandle> verts;
    rval = mbi_->get_connectivity(&tet, 1, verts);
    if (rval != moab::MB_SUCCESS) {
      fatal_error("Failed to get connectivity of tet on umesh: " + filename_);
    }

    moab::CartVect p[4];
    rval = mbi_->get_coords(verts.data(), verts.size(), p[0].array());
    if (rval != moab::MB_SUCCESS) {
      fatal_error("Failed to get coordinates of a tet in umesh: " + filename_);
    }

    moab::Matrix3 a(p[1] - p[0], p[2] - p[0], p[3] - p[0], true);

    // invert now to avoid this cost later
    a = a.transpose().inverse();
    baryc_data_.at(get_bin_from_ent_handle(tet)) = a;
  }
}

bool MOABMesh::point_in_tet(
  const moab::CartVect& r, moab::EntityHandle tet) const
{

  moab::ErrorCode rval;

  // get tet vertices
  vector<moab::EntityHandle> verts;
  rval = mbi_->get_connectivity(&tet, 1, verts);
  if (rval != moab::MB_SUCCESS) {
    warning("Failed to get vertices of tet in umesh: " + filename_);
    return false;
  }

  // first vertex is used as a reference point for the barycentric data -
  // retrieve its coordinates
  moab::CartVect p_zero;
  rval = mbi_->get_coords(verts.data(), 1, p_zero.array());
  if (rval != moab::MB_SUCCESS) {
    warning("Failed to get coordinates of a vertex in "
            "unstructured mesh: " +
            filename_);
    return false;
  }

  // look up barycentric data
  int idx = get_bin_from_ent_handle(tet);
  const moab::Matrix3& a_inv = baryc_data_[idx];

  moab::CartVect bary_coords = a_inv * (r - p_zero);

  return (bary_coords[0] >= 0.0 && bary_coords[1] >= 0.0 &&
          bary_coords[2] >= 0.0 &&
          bary_coords[0] + bary_coords[1] + bary_coords[2] <= 1.0);
}

int MOABMesh::get_bin_from_index(int idx) const
{
  if (idx >= n_bins()) {
    fatal_error(fmt::format("Invalid bin index: {}", idx));
  }
  return ehs_[idx] - ehs_[0];
}

int MOABMesh::get_index(const Position& r, bool* in_mesh) const
{
  int bin = get_bin(r);
  *in_mesh = bin != -1;
  return bin;
}

int MOABMesh::get_index_from_bin(int bin) const
{
  return bin;
}

std::pair<vector<double>, vector<double>> MOABMesh::plot(
  Position plot_ll, Position plot_ur) const
{
  // TODO: Implement mesh lines
  return {};
}

int MOABMesh::get_bin_from_ent_handle(moab::EntityHandle eh) const
{
  int bin = eh - ehs_[0];
  if (bin >= n_bins()) {
    fatal_error(fmt::format("Invalid bin: {}", bin));
  }
  return bin;
}

moab::EntityHandle MOABMesh::get_ent_handle_from_bin(int bin) const
{
  if (bin >= n_bins()) {
    fatal_error(fmt::format("Invalid bin index: ", bin));
  }
  return ehs_[0] + bin;
}

int MOABMesh::n_bins() const
{
  return ehs_.size();
}

int MOABMesh::n_surface_bins() const
{
  // collect all triangles in the set of tets for this mesh
  moab::Range tris;
  moab::ErrorCode rval;
  rval = mbi_->get_entities_by_type(0, moab::MBTRI, tris);
  if (rval != moab::MB_SUCCESS) {
    warning("Failed to get all triangles in the mesh instance");
    return -1;
  }
  return 2 * tris.size();
}

Position MOABMesh::centroid(int bin) const
{
  moab::ErrorCode rval;

  auto tet = this->get_ent_handle_from_bin(bin);

  // look up the tet connectivity
  vector<moab::EntityHandle> conn;
  rval = mbi_->get_connectivity(&tet, 1, conn);
  if (rval != moab::MB_SUCCESS) {
    warning("Failed to get connectivity of a mesh element.");
    return {};
  }

  // get the coordinates
  vector<moab::CartVect> coords(conn.size());
  rval = mbi_->get_coords(conn.data(), conn.size(), coords[0].array());
  if (rval != moab::MB_SUCCESS) {
    warning("Failed to get the coordinates of a mesh element.");
    return {};
  }

  // compute the centroid of the element vertices
  moab::CartVect centroid(0.0, 0.0, 0.0);
  for (const auto& coord : coords) {
    centroid += coord;
  }
  centroid /= double(coords.size());

  return {centroid[0], centroid[1], centroid[2]};
}

std::pair<moab::Tag, moab::Tag> MOABMesh::get_score_tags(
  std::string score) const
{
  moab::ErrorCode rval;
  // add a tag to the mesh
  // all scores are treated as a single value
  // with an uncertainty
  moab::Tag value_tag;

  // create the value tag if not present and get handle
  double default_val = 0.0;
  auto val_string = score + "_mean";
  rval = mbi_->tag_get_handle(val_string.c_str(), 1, moab::MB_TYPE_DOUBLE,
    value_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &default_val);
  if (rval != moab::MB_SUCCESS) {
    auto msg =
      fmt::format("Could not create or retrieve the value tag for the score {}"
                  " on unstructured mesh {}",
        score, id_);
    fatal_error(msg);
  }

  // create the std dev tag if not present and get handle
  moab::Tag error_tag;
  std::string err_string = score + "_std_dev";
  rval = mbi_->tag_get_handle(err_string.c_str(), 1, moab::MB_TYPE_DOUBLE,
    error_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &default_val);
  if (rval != moab::MB_SUCCESS) {
    auto msg =
      fmt::format("Could not create or retrieve the error tag for the score {}"
                  " on unstructured mesh {}",
        score, id_);
    fatal_error(msg);
  }

  // return the populated tag handles
  return {value_tag, error_tag};
}

void MOABMesh::add_score(const std::string& score)
{
  auto score_tags = get_score_tags(score);
  tag_names_.push_back(score);
}

void MOABMesh::remove_scores()
{
  for (const auto& name : tag_names_) {
    auto value_name = name + "_mean";
    moab::Tag tag;
    moab::ErrorCode rval = mbi_->tag_get_handle(value_name.c_str(), tag);
    if (rval != moab::MB_SUCCESS)
      return;

    rval = mbi_->tag_delete(tag);
    if (rval != moab::MB_SUCCESS) {
      auto msg = fmt::format("Failed to delete mesh tag for the score {}"
                             " on unstructured mesh {}",
        name, id_);
      fatal_error(msg);
    }

    auto std_dev_name = name + "_std_dev";
    rval = mbi_->tag_get_handle(std_dev_name.c_str(), tag);
    if (rval != moab::MB_SUCCESS) {
      auto msg =
        fmt::format("Std. Dev. mesh tag does not exist for the score {}"
                    " on unstructured mesh {}",
          name, id_);
    }

    rval = mbi_->tag_delete(tag);
    if (rval != moab::MB_SUCCESS) {
      auto msg = fmt::format("Failed to delete mesh tag for the score {}"
                             " on unstructured mesh {}",
        name, id_);
      fatal_error(msg);
    }
  }
  tag_names_.clear();
}

void MOABMesh::set_score_data(const std::string& score,
  const vector<double>& values, const vector<double>& std_dev)
{
  auto score_tags = this->get_score_tags(score);

  moab::ErrorCode rval;
  // set the score value
  rval = mbi_->tag_set_data(score_tags.first, ehs_, values.data());
  if (rval != moab::MB_SUCCESS) {
    auto msg = fmt::format("Failed to set the tally value for score '{}' "
                           "on unstructured mesh {}",
      score, id_);
    warning(msg);
  }

  // set the error value
  rval = mbi_->tag_set_data(score_tags.second, ehs_, std_dev.data());
  if (rval != moab::MB_SUCCESS) {
    auto msg = fmt::format("Failed to set the tally error for score '{}' "
                           "on unstructured mesh {}",
      score, id_);
    warning(msg);
  }
}

void MOABMesh::write(const std::string& base_filename) const
{
  // add extension to the base name
  auto filename = base_filename + ".vtk";
  write_message(5, "Writing unstructured mesh {}...", filename);
  filename = settings::path_output + filename;

  // write the tetrahedral elements of the mesh only
  // to avoid clutter from zero-value data on other
  // elements during visualization
  moab::ErrorCode rval;
  rval = mbi_->write_mesh(filename.c_str(), &tetset_, 1);
  if (rval != moab::MB_SUCCESS) {
    auto msg = fmt::format("Failed to write unstructured mesh {}", id_);
    warning(msg);
  }
}

#endif

#ifdef LIBMESH

LibMesh::LibMesh(pugi::xml_node node) : UnstructuredMesh(node)
{
  initialize();
}

LibMesh::LibMesh(const std::string& filename)
{
  filename_ = filename;
  initialize();
}

void LibMesh::initialize()
{
  if (!settings::libmesh_comm) {
    fatal_error(
      "Attempting to use an unstructured mesh without a libMesh communicator.");
  }

  // assuming that unstructured meshes used in OpenMC are 3D
  n_dimension_ = 3;

  m_ = make_unique<libMesh::Mesh>(*settings::libmesh_comm, n_dimension_);
  m_->read(filename_);
  m_->prepare_for_use();

  // ensure that the loaded mesh is 3 dimensional
  if (m_->mesh_dimension() != n_dimension_) {
    fatal_error(fmt::format("Mesh file {} specified for use in an unstructured "
                            "mesh is not a 3D mesh.",
      filename_));
  }

  // create an equation system for storing values
  eq_system_name_ = fmt::format("mesh_{}_system", id_);

  equation_systems_ = make_unique<libMesh::EquationSystems>(*m_);
  libMesh::ExplicitSystem& eq_sys =
    equation_systems_->add_system<libMesh::ExplicitSystem>(eq_system_name_);

#ifdef _OPENMP
  int n_threads = omp_get_max_threads();
#else
  int n_threads = 1;
#endif

  for (int i = 0; i < n_threads; i++) {
    pl_.emplace_back(m_->sub_point_locator());
    pl_.back()->set_contains_point_tol(FP_COINCIDENT);
    pl_.back()->enable_out_of_mesh_mode();
  }

  // store first element in the mesh to use as an offset for bin indices
  auto first_elem = *m_->elements_begin();
  first_element_id_ = first_elem->id();

  // bounding box for the mesh for quick rejection checks
  bbox_ = libMesh::MeshTools::create_bounding_box(*m_);
}

Position LibMesh::centroid(int bin) const
{
  const auto& elem = this->get_element_from_bin(bin);
  auto centroid = elem.centroid();
  return {centroid(0), centroid(1), centroid(2)};
}

std::string LibMesh::library() const
{
  return "libmesh";
}

int LibMesh::n_bins() const
{
  return m_->n_elem();
}

int LibMesh::n_surface_bins() const
{
  int n_bins = 0;
  for (int i = 0; i < this->n_bins(); i++) {
    const libMesh::Elem& e = get_element_from_bin(i);
    n_bins += e.n_faces();
    // if this is a boundary element, it will only be visited once,
    // the number of surface bins is incremented to
    for (auto neighbor_ptr : e.neighbor_ptr_range()) {
      // null neighbor pointer indicates a boundary face
      if (!neighbor_ptr) {
        n_bins++;
      }
    }
  }
  return n_bins;
}

void LibMesh::add_score(const std::string& var_name)
{
  // check if this is a new variable
  std::string value_name = var_name + "_mean";
  if (!variable_map_.count(value_name)) {
    auto& eqn_sys = equation_systems_->get_system(eq_system_name_);
    auto var_num =
      eqn_sys.add_variable(value_name, libMesh::CONSTANT, libMesh::MONOMIAL);
    variable_map_[value_name] = var_num;
  }

  std::string std_dev_name = var_name + "_std_dev";
  // check if this is a new variable
  if (!variable_map_.count(std_dev_name)) {
    auto& eqn_sys = equation_systems_->get_system(eq_system_name_);
    auto var_num =
      eqn_sys.add_variable(std_dev_name, libMesh::CONSTANT, libMesh::MONOMIAL);
    variable_map_[std_dev_name] = var_num;
  }
}

void LibMesh::remove_scores()
{
  auto& eqn_sys = equation_systems_->get_system(eq_system_name_);
  eqn_sys.clear();
  variable_map_.clear();
}

void LibMesh::set_score_data(const std::string& var_name,
  const vector<double>& values, const vector<double>& std_dev)
{
  auto& eqn_sys = equation_systems_->get_system(eq_system_name_);

  if (!eqn_sys.is_initialized()) {
    equation_systems_->init();
  }

  const libMesh::DofMap& dof_map = eqn_sys.get_dof_map();

  // look up the value variable
  std::string value_name = var_name + "_mean";
  unsigned int value_num = variable_map_.at(value_name);
  // look up the std dev variable
  std::string std_dev_name = var_name + "_std_dev";
  unsigned int std_dev_num = variable_map_.at(std_dev_name);

  for (auto it = m_->local_elements_begin(); it != m_->local_elements_end();
       it++) {
    auto bin = get_bin_from_element(*it);

    // set value
    vector<libMesh::dof_id_type> value_dof_indices;
    dof_map.dof_indices(*it, value_dof_indices, value_num);
    Ensures(value_dof_indices.size() == 1);
    eqn_sys.solution->set(value_dof_indices[0], values.at(bin));

    // set std dev
    vector<libMesh::dof_id_type> std_dev_dof_indices;
    dof_map.dof_indices(*it, std_dev_dof_indices, std_dev_num);
    Ensures(std_dev_dof_indices.size() == 1);
    eqn_sys.solution->set(std_dev_dof_indices[0], std_dev.at(bin));
  }
}

void LibMesh::write(const std::string& filename) const
{
  write_message(fmt::format(
    "Writing file: {}.e for unstructured mesh {}", filename, this->id_));
  libMesh::ExodusII_IO exo(*m_);
  std::set<std::string> systems_out = {eq_system_name_};
  exo.write_discontinuous_exodusII(
    filename + ".e", *equation_systems_, &systems_out);
}

void LibMesh::bins_crossed(Position r0, Position r1, const Direction& u,
  vector<int>& bins, vector<double>& lengths) const
{
  // TODO: Implement triangle crossings here
  fatal_error("Tracklength tallies on libMesh instances are not implemented.");
}

int LibMesh::get_bin(Position r) const
{
  // look-up a tet using the point locator
  libMesh::Point p(r.x, r.y, r.z);

  // quick rejection check
  if (!bbox_.contains_point(p)) {
    return -1;
  }

#ifdef _OPENMP
  int thread_num = omp_get_thread_num();
#else
  int thread_num = 0;
#endif

  const auto& point_locator = pl_.at(thread_num);

  const auto elem_ptr = (*point_locator)(p);
  return elem_ptr ? get_bin_from_element(elem_ptr) : -1;
}

int LibMesh::get_bin_from_element(const libMesh::Elem* elem) const
{
  int bin = elem->id() - first_element_id_;
  if (bin >= n_bins() || bin < 0) {
    fatal_error(fmt::format("Invalid bin: {}", bin));
  }
  return bin;
}

std::pair<vector<double>, vector<double>> LibMesh::plot(
  Position plot_ll, Position plot_ur) const
{
  return {};
}

const libMesh::Elem& LibMesh::get_element_from_bin(int bin) const
{
  return m_->elem_ref(bin);
}

double LibMesh::volume(int bin) const
{
  return m_->elem_ref(bin).volume();
}

#endif // LIBMESH

//==============================================================================
// Non-member functions
//==============================================================================

void read_meshes(pugi::xml_node root)
{
  for (auto node : root.children("mesh")) {
    std::string mesh_type;
    if (check_for_node(node, "type")) {
      mesh_type = get_node_value(node, "type", true, true);
    } else {
      mesh_type = "regular";
    }

    // determine the mesh library to use
    std::string mesh_lib;
    if (check_for_node(node, "library")) {
      mesh_lib = get_node_value(node, "library", true, true);
    }

    // Read mesh and add to vector
    if (mesh_type == "regular") {
      model::meshes.push_back(make_unique<RegularMesh>(node));
    } else if (mesh_type == "rectilinear") {
      model::meshes.push_back(make_unique<RectilinearMesh>(node));
#ifdef DAGMC
    } else if (mesh_type == "unstructured" && mesh_lib == "moab") {
      model::meshes.push_back(make_unique<MOABMesh>(node));
#endif
#ifdef LIBMESH
    } else if (mesh_type == "unstructured" && mesh_lib == "libmesh") {
      model::meshes.push_back(make_unique<LibMesh>(node));
#endif
    } else if (mesh_type == "unstructured") {
      fatal_error("Unstructured mesh support is not enabled or the mesh "
                  "library is invalid.");
    } else {
      fatal_error("Invalid mesh type: " + mesh_type);
    }

    // Map ID to position in vector
    model::mesh_map[model::meshes.back()->id_] = model::meshes.size() - 1;
  }
}

void meshes_to_hdf5(hid_t group)
{
  // Write number of meshes
  hid_t meshes_group = create_group(group, "meshes");
  int32_t n_meshes = model::meshes.size();
  write_attribute(meshes_group, "n_meshes", n_meshes);

  if (n_meshes > 0) {
    // Write IDs of meshes
    vector<int> ids;
    for (const auto& m : model::meshes) {
      m->to_hdf5(meshes_group);
      ids.push_back(m->id_);
    }
    write_attribute(meshes_group, "ids", ids);
  }

  close_group(meshes_group);
}

void free_memory_mesh()
{
  model::meshes.clear();
  model::mesh_map.clear();
}

extern "C" int n_meshes()
{
  return model::meshes.size();
}

} // namespace openmc
