#include "openmc/mesh.h"
#include <algorithm> // for copy, equal, min, min_element
#include <cmath>     // for ceil
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

#ifdef LIBMESH
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/numeric_vector.h"
#endif

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
  // Read mesh id
  id_ = std::stoi(get_node_value(node, "id"));
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

vector<double> Mesh::volumes() const
{
  vector<double> volumes(n_bins());
  for (int i = 0; i < n_bins(); i++) {
    volumes[i] = this->volume(i);
  }
  return volumes;
}

int Mesh::material_volumes(
  int n_sample, int bin, gsl::span<MaterialVolume> result, uint64_t* seed) const
{
  vector<int32_t> materials;
  vector<int64_t> hits;

#pragma omp parallel
  {
    vector<int32_t> local_materials;
    vector<int64_t> local_hits;
    GeometryState geom;

#pragma omp for
    for (int i = 0; i < n_sample; ++i) {
      // Get seed for i-th sample
      uint64_t seed_i = future_seed(3 * i, *seed);

      // Sample position and set geometry state
      geom.r() = this->sample_element(bin, &seed_i);
      geom.u() = {1., 0., 0.};
      geom.n_coord() = 1;

      // If this location is not in the geometry at all, move on to next block
      if (!exhaustive_find_cell(geom))
        continue;

      int i_material = geom.material();

      // Check if this material was previously hit and if so, increment count
      auto it =
        std::find(local_materials.begin(), local_materials.end(), i_material);
      if (it == local_materials.end()) {
        local_materials.push_back(i_material);
        local_hits.push_back(1);
      } else {
        local_hits[it - local_materials.begin()]++;
      }
    } // omp for

    // Reduce index/hits lists from each thread into a single copy
    reduce_indices_hits(local_materials, local_hits, materials, hits);
  } // omp parallel

  // Advance RNG seed
  advance_prn_seed(3 * n_sample, seed);

  // Make sure span passed in is large enough
  if (hits.size() > result.size()) {
    return -1;
  }

  // Convert hits to fractions
  for (int i_mat = 0; i_mat < hits.size(); ++i_mat) {
    double fraction = double(hits[i_mat]) / n_sample;
    result[i_mat].material = materials[i_mat];
    result[i_mat].volume = fraction * this->volume(bin);
  }
  return hits.size();
}

vector<Mesh::MaterialVolume> Mesh::material_volumes(
  int n_sample, int bin, uint64_t* seed) const
{
  // Create result vector with space for 8 pairs
  vector<Mesh::MaterialVolume> result;
  result.reserve(8);

  int size = -1;
  while (true) {
    // Get material volumes
    size = this->material_volumes(
      n_sample, bin, {result.data(), result.data() + result.capacity()}, seed);

    // If capacity was sufficient, resize the vector and return
    if (size >= 0) {
      result.resize(size);
      break;
    }

    // Otherwise, increase capacity of the vector
    result.reserve(2 * result.capacity());
  }

  return result;
}

//==============================================================================
// Structured Mesh implementation
//==============================================================================

std::string StructuredMesh::bin_label(int bin) const
{
  MeshIndex ijk = get_indices_from_bin(bin);

  if (n_dimension_ > 2) {
    return fmt::format("Mesh Index ({}, {}, {})", ijk[0], ijk[1], ijk[2]);
  } else if (n_dimension_ > 1) {
    return fmt::format("Mesh Index ({}, {})", ijk[0], ijk[1]);
  } else {
    return fmt::format("Mesh Index ({})", ijk[0]);
  }
}

xt::xtensor<int, 1> StructuredMesh::get_x_shape() const
{
  // because method is const, shape_ is const as well and can't be adapted
  auto tmp_shape = shape_;
  return xt::adapt(tmp_shape, {n_dimension_});
}

Position StructuredMesh::sample_element(
  const MeshIndex& ijk, uint64_t* seed) const
{
  // lookup the lower/upper bounds for the mesh element
  double x_min = negative_grid_boundary(ijk, 0);
  double x_max = positive_grid_boundary(ijk, 0);

  double y_min = (n_dimension_ >= 2) ? negative_grid_boundary(ijk, 1) : 0.0;
  double y_max = (n_dimension_ >= 2) ? positive_grid_boundary(ijk, 1) : 0.0;

  double z_min = (n_dimension_ == 3) ? negative_grid_boundary(ijk, 2) : 0.0;
  double z_max = (n_dimension_ == 3) ? positive_grid_boundary(ijk, 2) : 0.0;

  return {x_min + (x_max - x_min) * prn(seed),
    y_min + (y_max - y_min) * prn(seed), z_min + (z_max - z_min) * prn(seed)};
}

//==============================================================================
// Unstructured Mesh implementation
//==============================================================================

UnstructuredMesh::UnstructuredMesh(pugi::xml_node node) : Mesh(node)
{

  // check the mesh type
  if (check_for_node(node, "type")) {
    auto temp = get_node_value(node, "type", true, true);
    if (temp != mesh_type) {
      fatal_error(fmt::format("Invalid mesh type: {}", temp));
    }
  }

  // check if a length unit multiplier was specified
  if (check_for_node(node, "length_multiplier")) {
    length_multiplier_ = std::stod(get_node_value(node, "length_multiplier"));
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

  if (check_for_node(node, "options")) {
    options_ = get_node_value(node, "options");
  }

  // check if mesh tally data should be written with
  // statepoint files
  if (check_for_node(node, "output")) {
    output_ = get_node_value_bool(node, "output");
  }
}

void UnstructuredMesh::determine_bounds()
{
  double xmin = INFTY;
  double ymin = INFTY;
  double zmin = INFTY;
  double xmax = -INFTY;
  double ymax = -INFTY;
  double zmax = -INFTY;
  int n = this->n_vertices();
  for (int i = 0; i < n; ++i) {
    auto v = this->vertex(i);
    xmin = std::min(v.x, xmin);
    ymin = std::min(v.y, ymin);
    zmin = std::min(v.z, zmin);
    xmax = std::max(v.x, xmax);
    ymax = std::max(v.y, ymax);
    zmax = std::max(v.z, zmax);
  }
  lower_left_ = {xmin, ymin, zmin};
  upper_right_ = {xmax, ymax, zmax};
}

Position UnstructuredMesh::sample_tet(
  std::array<Position, 4> coords, uint64_t* seed) const
{
  // Uniform distribution
  double s = prn(seed);
  double t = prn(seed);
  double u = prn(seed);

  // From PyNE implementation of moab tet sampling C. Rocchini & P. Cignoni
  // (2000) Generating Random Points in a Tetrahedron, Journal of Graphics
  // Tools, 5:4, 9-12, DOI: 10.1080/10867651.2000.10487528
  if (s + t > 1) {
    s = 1.0 - s;
    t = 1.0 - t;
  }
  if (s + t + u > 1) {
    if (t + u > 1) {
      double old_t = t;
      t = 1.0 - u;
      u = 1.0 - s - old_t;
    } else if (t + u <= 1) {
      double old_s = s;
      s = 1.0 - t - u;
      u = old_s + t + u - 1;
    }
  }
  return s * (coords[1] - coords[0]) + t * (coords[2] - coords[0]) +
         u * (coords[3] - coords[0]) + coords[0];
}

const std::string UnstructuredMesh::mesh_type = "unstructured";

std::string UnstructuredMesh::get_mesh_type() const
{
  return mesh_type;
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

  write_dataset(mesh_group, "type", mesh_type);
  write_dataset(mesh_group, "filename", filename_);
  write_dataset(mesh_group, "library", this->library());
  if (!options_.empty()) {
    write_attribute(mesh_group, "options", options_);
  }

  if (length_multiplier_ > 0.0)
    write_dataset(mesh_group, "length_multiplier", length_multiplier_);

  // write vertex coordinates
  xt::xtensor<double, 2> vertices({static_cast<size_t>(this->n_vertices()), 3});
  for (int i = 0; i < this->n_vertices(); i++) {
    auto v = this->vertex(i);
    xt::view(vertices, i, xt::all()) = xt::xarray<double>({v.x, v.y, v.z});
  }
  write_dataset(mesh_group, "vertices", vertices);

  int num_elem_skipped = 0;

  // write element types and connectivity
  vector<double> volumes;
  xt::xtensor<int, 2> connectivity({static_cast<size_t>(this->n_bins()), 8});
  xt::xtensor<int, 2> elem_types({static_cast<size_t>(this->n_bins()), 1});
  for (int i = 0; i < this->n_bins(); i++) {
    auto conn = this->connectivity(i);

    volumes.emplace_back(this->volume(i));

    // write linear tet element
    if (conn.size() == 4) {
      xt::view(elem_types, i, xt::all()) =
        static_cast<int>(ElementType::LINEAR_TET);
      xt::view(connectivity, i, xt::all()) =
        xt::xarray<int>({conn[0], conn[1], conn[2], conn[3], -1, -1, -1, -1});
      // write linear hex element
    } else if (conn.size() == 8) {
      xt::view(elem_types, i, xt::all()) =
        static_cast<int>(ElementType::LINEAR_HEX);
      xt::view(connectivity, i, xt::all()) = xt::xarray<int>({conn[0], conn[1],
        conn[2], conn[3], conn[4], conn[5], conn[6], conn[7]});
    } else {
      num_elem_skipped++;
      xt::view(elem_types, i, xt::all()) =
        static_cast<int>(ElementType::UNSUPPORTED);
      xt::view(connectivity, i, xt::all()) = -1;
    }
  }

  // warn users that some elements were skipped
  if (num_elem_skipped > 0) {
    warning(fmt::format("The connectivity of {} elements "
                        "on mesh {} were not written "
                        "because they are not of type linear tet/hex.",
      num_elem_skipped, this->id_));
  }

  write_dataset(mesh_group, "volumes", volumes);
  write_dataset(mesh_group, "connectivity", connectivity);
  write_dataset(mesh_group, "element_types", elem_types);

  close_group(mesh_group);
}

void UnstructuredMesh::set_length_multiplier(double length_multiplier)
{
  length_multiplier_ = length_multiplier;
}

ElementType UnstructuredMesh::element_type(int bin) const
{
  auto conn = connectivity(bin);

  if (conn.size() == 4)
    return ElementType::LINEAR_TET;
  else if (conn.size() == 8)
    return ElementType::LINEAR_HEX;
  else
    return ElementType::UNSUPPORTED;
}

StructuredMesh::MeshIndex StructuredMesh::get_indices(
  Position r, bool& in_mesh) const
{
  MeshIndex ijk;
  in_mesh = true;
  for (int i = 0; i < n_dimension_; ++i) {
    ijk[i] = get_index_in_direction(r[i], i);

    if (ijk[i] < 1 || ijk[i] > shape_[i])
      in_mesh = false;
  }
  return ijk;
}

int StructuredMesh::get_bin_from_indices(const MeshIndex& ijk) const
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

StructuredMesh::MeshIndex StructuredMesh::get_indices_from_bin(int bin) const
{
  MeshIndex ijk;
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
  return ijk;
}

int StructuredMesh::get_bin(Position r) const
{
  // Determine indices
  bool in_mesh;
  MeshIndex ijk = get_indices(r, in_mesh);
  if (!in_mesh)
    return -1;

  // Convert indices to bin
  return get_bin_from_indices(ijk);
}

int StructuredMesh::n_bins() const
{
  return std::accumulate(
    shape_.begin(), shape_.begin() + n_dimension_, 1, std::multiplies<>());
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

// raytrace through the mesh. The template class T will do the tallying.
// A modern optimizing compiler can recognize the noop method of T and eleminate
// that call entirely.
template<class T>
void StructuredMesh::raytrace_mesh(
  Position r0, Position r1, const Direction& u, T tally) const
{
  // TODO: when c++-17 is available, use "if constexpr ()" to compile-time
  // enable/disable tally calls for now, T template type needs to provide both
  // surface and track methods, which might be empty. modern optimizing
  // compilers will (hopefully) eliminate the complete code (including
  // calculation of parameters) but for the future: be explicit

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

void StructuredMesh::bins_crossed(Position r0, Position r1, const Direction& u,
  vector<int>& bins, vector<double>& lengths) const
{

  // Helper tally class.
  // stores a pointer to the mesh class and references to bins and lengths
  // parameters. Performs the actual tally through the track method.
  struct TrackAggregator {
    TrackAggregator(
      const StructuredMesh* _mesh, vector<int>& _bins, vector<double>& _lengths)
      : mesh(_mesh), bins(_bins), lengths(_lengths)
    {}
    void surface(const MeshIndex& ijk, int k, bool max, bool inward) const {}
    void track(const MeshIndex& ijk, double l) const
    {
      bins.push_back(mesh->get_bin_from_indices(ijk));
      lengths.push_back(l);
    }

    const StructuredMesh* mesh;
    vector<int>& bins;
    vector<double>& lengths;
  };

  // Perform the mesh raytrace with the helper class.
  raytrace_mesh(r0, r1, u, TrackAggregator(this, bins, lengths));
}

void StructuredMesh::surface_bins_crossed(
  Position r0, Position r1, const Direction& u, vector<int>& bins) const
{

  // Helper tally class.
  // stores a pointer to the mesh class and a reference to the bins parameter.
  // Performs the actual tally through the surface method.
  struct SurfaceAggregator {
    SurfaceAggregator(const StructuredMesh* _mesh, vector<int>& _bins)
      : mesh(_mesh), bins(_bins)
    {}
    void surface(const MeshIndex& ijk, int k, bool max, bool inward) const
    {
      int i_bin =
        4 * mesh->n_dimension_ * mesh->get_bin_from_indices(ijk) + 4 * k;
      if (max)
        i_bin += 2;
      if (inward)
        i_bin += 1;
      bins.push_back(i_bin);
    }
    void track(const MeshIndex& idx, double l) const {}

    const StructuredMesh* mesh;
    vector<int>& bins;
  };

  // Perform the mesh raytrace with the helper class.
  raytrace_mesh(r0, r1, u, SurfaceAggregator(this, bins));
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

  xt::xtensor<int, 1> shape = get_node_xarray<int>(node, "dimension");
  int n = n_dimension_ = shape.size();
  if (n != 1 && n != 2 && n != 3) {
    fatal_error("Mesh must be one, two, or three dimensions.");
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
    fatal_error("Must specify <lower_left> on a mesh.");
  }

  // Make sure lower_left and dimension match
  if (shape.size() != lower_left_.size()) {
    fatal_error("Number of entries on <lower_left> must be the same "
                "as the number of entries on <dimension>.");
  }

  if (check_for_node(node, "width")) {
    // Make sure one of upper-right or width were specified
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
                  "the <lower_left> coordinates on a tally mesh.");
    }

    // Set width
    width_ = xt::eval((upper_right_ - lower_left_) / shape);
  } else {
    fatal_error("Must specify either <upper_right> or <width> on a mesh.");
  }

  // Set material volumes
  volume_frac_ = 1.0 / xt::prod(shape)();

  element_volume_ = 1.0;
  for (int i = 0; i < n_dimension_; i++) {
    element_volume_ *= width_[i];
  }
}

int RegularMesh::get_index_in_direction(double r, int i) const
{
  return std::ceil((r - lower_left_[i]) / width_[i]);
}

const std::string RegularMesh::mesh_type = "regular";

std::string RegularMesh::get_mesh_type() const
{
  return mesh_type;
}

double RegularMesh::positive_grid_boundary(const MeshIndex& ijk, int i) const
{
  return lower_left_[i] + ijk[i] * width_[i];
}

double RegularMesh::negative_grid_boundary(const MeshIndex& ijk, int i) const
{
  return lower_left_[i] + (ijk[i] - 1) * width_[i];
}

StructuredMesh::MeshDistance RegularMesh::distance_to_grid_boundary(
  const MeshIndex& ijk, int i, const Position& r0, const Direction& u,
  double l) const
{
  MeshDistance d;
  d.next_index = ijk[i];
  if (std::abs(u[i]) < FP_PRECISION)
    return d;

  d.max_surface = (u[i] > 0);
  if (d.max_surface && (ijk[i] <= shape_[i])) {
    d.next_index++;
    d.distance = (positive_grid_boundary(ijk, i) - r0[i]) / u[i];
  } else if (!d.max_surface && (ijk[i] >= 1)) {
    d.next_index--;
    d.distance = (negative_grid_boundary(ijk, i) - r0[i]) / u[i];
  }
  return d;
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
  write_dataset(mesh_group, "dimension", get_x_shape());
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

double RegularMesh::volume(const MeshIndex& ijk) const
{
  return element_volume_;
}

//==============================================================================
// RectilinearMesh implementation
//==============================================================================

RectilinearMesh::RectilinearMesh(pugi::xml_node node) : StructuredMesh {node}
{
  n_dimension_ = 3;

  grid_[0] = get_node_array<double>(node, "x_grid");
  grid_[1] = get_node_array<double>(node, "y_grid");
  grid_[2] = get_node_array<double>(node, "z_grid");

  if (int err = set_grid()) {
    fatal_error(openmc_err_msg);
  }
}

const std::string RectilinearMesh::mesh_type = "rectilinear";

std::string RectilinearMesh::get_mesh_type() const
{
  return mesh_type;
}

double RectilinearMesh::positive_grid_boundary(
  const MeshIndex& ijk, int i) const
{
  return grid_[i][ijk[i]];
}

double RectilinearMesh::negative_grid_boundary(
  const MeshIndex& ijk, int i) const
{
  return grid_[i][ijk[i] - 1];
}

StructuredMesh::MeshDistance RectilinearMesh::distance_to_grid_boundary(
  const MeshIndex& ijk, int i, const Position& r0, const Direction& u,
  double l) const
{
  MeshDistance d;
  d.next_index = ijk[i];
  if (std::abs(u[i]) < FP_PRECISION)
    return d;

  d.max_surface = (u[i] > 0);
  if (d.max_surface && (ijk[i] <= shape_[i])) {
    d.next_index++;
    d.distance = (positive_grid_boundary(ijk, i) - r0[i]) / u[i];
  } else if (!d.max_surface && (ijk[i] > 0)) {
    d.next_index--;
    d.distance = (negative_grid_boundary(ijk, i) - r0[i]) / u[i];
  }
  return d;
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
    if (std::adjacent_find(g.begin(), g.end(), std::greater_equal<>()) !=
        g.end()) {
      set_errmsg("Values in for x-, y-, and z- grids for "
                 "rectilinear meshes must be sorted and unique.");
      return OPENMC_E_INVALID_ARGUMENT;
    }
  }

  lower_left_ = {grid_[0].front(), grid_[1].front(), grid_[2].front()};
  upper_right_ = {grid_[0].back(), grid_[1].back(), grid_[2].back()};

  return 0;
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

double RectilinearMesh::volume(const MeshIndex& ijk) const
{
  double vol {1.0};

  for (int i = 0; i < n_dimension_; i++) {
    vol *= grid_[i][ijk[i]] - grid_[i][ijk[i] - 1];
  }
  return vol;
}

//==============================================================================
// CylindricalMesh implementation
//==============================================================================

CylindricalMesh::CylindricalMesh(pugi::xml_node node)
  : PeriodicStructuredMesh {node}
{
  n_dimension_ = 3;
  grid_[0] = get_node_array<double>(node, "r_grid");
  grid_[1] = get_node_array<double>(node, "phi_grid");
  grid_[2] = get_node_array<double>(node, "z_grid");
  origin_ = get_node_position(node, "origin");

  if (int err = set_grid()) {
    fatal_error(openmc_err_msg);
  }
}

const std::string CylindricalMesh::mesh_type = "cylindrical";

std::string CylindricalMesh::get_mesh_type() const
{
  return mesh_type;
}

StructuredMesh::MeshIndex CylindricalMesh::get_indices(
  Position r, bool& in_mesh) const
{
  local_coords(r);

  Position mapped_r;
  mapped_r[0] = std::hypot(r.x, r.y);
  mapped_r[2] = r[2];

  if (mapped_r[0] < FP_PRECISION) {
    mapped_r[1] = 0.0;
  } else {
    mapped_r[1] = std::atan2(r.y, r.x);
    if (mapped_r[1] < 0)
      mapped_r[1] += 2 * M_PI;
  }

  MeshIndex idx = StructuredMesh::get_indices(mapped_r, in_mesh);

  idx[1] = sanitize_phi(idx[1]);

  return idx;
}

Position CylindricalMesh::sample_element(
  const MeshIndex& ijk, uint64_t* seed) const
{
  double r_min = this->r(ijk[0] - 1);
  double r_max = this->r(ijk[0]);

  double phi_min = this->phi(ijk[1] - 1);
  double phi_max = this->phi(ijk[1]);

  double z_min = this->z(ijk[2] - 1);
  double z_max = this->z(ijk[2]);

  double r_min_sq = r_min * r_min;
  double r_max_sq = r_max * r_max;
  double r = std::sqrt(uniform_distribution(r_min_sq, r_max_sq, seed));
  double phi = uniform_distribution(phi_min, phi_max, seed);
  double z = uniform_distribution(z_min, z_max, seed);

  double x = r * std::cos(phi);
  double y = r * std::sin(phi);

  return origin_ + Position(x, y, z);
}

double CylindricalMesh::find_r_crossing(
  const Position& r, const Direction& u, double l, int shell) const
{

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

double CylindricalMesh::find_phi_crossing(
  const Position& r, const Direction& u, double l, int shell) const
{
  // Phi grid is [0, 2π], thus there is no real surface to cross
  if (full_phi_ && (shape_[1] == 1))
    return INFTY;

  shell = sanitize_phi(shell);

  const double p0 = grid_[1][shell];

  // solve y(s)/x(s) = tan(p0) = sin(p0)/cos(p0)
  // => x(s) * cos(p0) = y(s) * sin(p0)
  // => (y + s * v) * cos(p0) = (x + s * u) * sin(p0)
  // = s * (v * cos(p0) - u * sin(p0)) = - (y * cos(p0) - x * sin(p0))

  const double c0 = std::cos(p0);
  const double s0 = std::sin(p0);

  const double denominator = (u.x * s0 - u.y * c0);

  // Check if direction of flight is not parallel to phi surface
  if (std::abs(denominator) > FP_PRECISION) {
    const double s = -(r.x * s0 - r.y * c0) / denominator;
    // Check if solution is in positive direction of flight and crosses the
    // correct phi surface (not -phi)
    if ((s > l) && ((c0 * (r.x + s * u.x) + s0 * (r.y + s * u.y)) > 0.0))
      return s;
  }

  return INFTY;
}

StructuredMesh::MeshDistance CylindricalMesh::find_z_crossing(
  const Position& r, const Direction& u, double l, int shell) const
{
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

StructuredMesh::MeshDistance CylindricalMesh::distance_to_grid_boundary(
  const MeshIndex& ijk, int i, const Position& r0, const Direction& u,
  double l) const
{
  Position r = r0 - origin_;

  if (i == 0) {

    return std::min(
      MeshDistance(ijk[i] + 1, true, find_r_crossing(r, u, l, ijk[i])),
      MeshDistance(ijk[i] - 1, false, find_r_crossing(r, u, l, ijk[i] - 1)));

  } else if (i == 1) {

    return std::min(MeshDistance(sanitize_phi(ijk[i] + 1), true,
                      find_phi_crossing(r, u, l, ijk[i])),
      MeshDistance(sanitize_phi(ijk[i] - 1), false,
        find_phi_crossing(r, u, l, ijk[i] - 1)));

  } else {
    return find_z_crossing(r, u, l, ijk[i]);
  }
}

int CylindricalMesh::set_grid()
{
  shape_ = {static_cast<int>(grid_[0].size()) - 1,
    static_cast<int>(grid_[1].size()) - 1,
    static_cast<int>(grid_[2].size()) - 1};

  for (const auto& g : grid_) {
    if (g.size() < 2) {
      set_errmsg("r-, phi-, and z- grids for cylindrical meshes "
                 "must each have at least 2 points");
      return OPENMC_E_INVALID_ARGUMENT;
    }
    if (std::adjacent_find(g.begin(), g.end(), std::greater_equal<>()) !=
        g.end()) {
      set_errmsg("Values in for r-, phi-, and z- grids for "
                 "cylindrical meshes must be sorted and unique.");
      return OPENMC_E_INVALID_ARGUMENT;
    }
  }
  if (grid_[0].front() < 0.0) {
    set_errmsg("r-grid for "
               "cylindrical meshes must start at r >= 0.");
    return OPENMC_E_INVALID_ARGUMENT;
  }
  if (grid_[1].front() < 0.0) {
    set_errmsg("phi-grid for "
               "cylindrical meshes must start at phi >= 0.");
    return OPENMC_E_INVALID_ARGUMENT;
  }
  if (grid_[1].back() > 2.0 * PI) {
    set_errmsg("phi-grids for "
               "cylindrical meshes must end with theta <= 2*pi.");

    return OPENMC_E_INVALID_ARGUMENT;
  }

  full_phi_ = (grid_[1].front() == 0.0) && (grid_[1].back() == 2.0 * PI);

  lower_left_ = {origin_[0] - grid_[0].back(), origin_[1] - grid_[0].back(),
    origin_[2] + grid_[2].front()};
  upper_right_ = {origin_[0] + grid_[0].back(), origin_[1] + grid_[0].back(),
    origin_[2] + grid_[2].back()};

  return 0;
}

int CylindricalMesh::get_index_in_direction(double r, int i) const
{
  return lower_bound_index(grid_[i].begin(), grid_[i].end(), r) + 1;
}

std::pair<vector<double>, vector<double>> CylindricalMesh::plot(
  Position plot_ll, Position plot_ur) const
{
  fatal_error("Plot of cylindrical Mesh not implemented");

  // Figure out which axes lie in the plane of the plot.
  array<vector<double>, 2> axis_lines;
  return {axis_lines[0], axis_lines[1]};
}

void CylindricalMesh::to_hdf5(hid_t group) const
{
  hid_t mesh_group = create_group(group, "mesh " + std::to_string(id_));

  write_dataset(mesh_group, "type", "cylindrical");
  write_dataset(mesh_group, "r_grid", grid_[0]);
  write_dataset(mesh_group, "phi_grid", grid_[1]);
  write_dataset(mesh_group, "z_grid", grid_[2]);
  write_dataset(mesh_group, "origin", origin_);

  close_group(mesh_group);
}

double CylindricalMesh::volume(const MeshIndex& ijk) const
{
  double r_i = grid_[0][ijk[0] - 1];
  double r_o = grid_[0][ijk[0]];

  double phi_i = grid_[1][ijk[1] - 1];
  double phi_o = grid_[1][ijk[1]];

  double z_i = grid_[2][ijk[2] - 1];
  double z_o = grid_[2][ijk[2]];

  return 0.5 * (r_o * r_o - r_i * r_i) * (phi_o - phi_i) * (z_o - z_i);
}

//==============================================================================
// SphericalMesh implementation
//==============================================================================

SphericalMesh::SphericalMesh(pugi::xml_node node)
  : PeriodicStructuredMesh {node}
{
  n_dimension_ = 3;

  grid_[0] = get_node_array<double>(node, "r_grid");
  grid_[1] = get_node_array<double>(node, "theta_grid");
  grid_[2] = get_node_array<double>(node, "phi_grid");
  origin_ = get_node_position(node, "origin");

  if (int err = set_grid()) {
    fatal_error(openmc_err_msg);
  }
}

const std::string SphericalMesh::mesh_type = "spherical";

std::string SphericalMesh::get_mesh_type() const
{
  return mesh_type;
}

StructuredMesh::MeshIndex SphericalMesh::get_indices(
  Position r, bool& in_mesh) const
{
  local_coords(r);

  Position mapped_r;
  mapped_r[0] = r.norm();

  if (mapped_r[0] < FP_PRECISION) {
    mapped_r[1] = 0.0;
    mapped_r[2] = 0.0;
  } else {
    mapped_r[1] = std::acos(r.z / mapped_r.x);
    mapped_r[2] = std::atan2(r.y, r.x);
    if (mapped_r[2] < 0)
      mapped_r[2] += 2 * M_PI;
  }

  MeshIndex idx = StructuredMesh::get_indices(mapped_r, in_mesh);

  idx[1] = sanitize_theta(idx[1]);
  idx[2] = sanitize_phi(idx[2]);

  return idx;
}

Position SphericalMesh::sample_element(
  const MeshIndex& ijk, uint64_t* seed) const
{
  double r_min = this->r(ijk[0] - 1);
  double r_max = this->r(ijk[0]);

  double theta_min = this->theta(ijk[1] - 1);
  double theta_max = this->theta(ijk[1]);

  double phi_min = this->phi(ijk[2] - 1);
  double phi_max = this->phi(ijk[2]);

  double cos_theta = uniform_distribution(theta_min, theta_max, seed);
  double sin_theta = std::sin(std::acos(cos_theta));
  double phi = uniform_distribution(phi_min, phi_max, seed);
  double r_min_cub = std::pow(r_min, 3);
  double r_max_cub = std::pow(r_max, 3);
  // might be faster to do rejection here?
  double r = std::cbrt(uniform_distribution(r_min_cub, r_max_cub, seed));

  double x = r * std::cos(phi) * sin_theta;
  double y = r * std::sin(phi) * sin_theta;
  double z = r * cos_theta;

  return origin_ + Position(x, y, z);
}

double SphericalMesh::find_r_crossing(
  const Position& r, const Direction& u, double l, int shell) const
{
  if ((shell < 0) || (shell > shape_[0]))
    return INFTY;

  // solve |r+s*u| = r0
  // |r+s*u| = |r| + 2*s*r*u + s^2 (|u|==1 !)
  const double r0 = grid_[0][shell];
  if (r0 == 0.0)
    return INFTY;
  const double p = r.dot(u);
  double c = r.dot(r) - r0 * r0;
  double D = p * p - c;

  if (std::abs(c) <= RADIAL_MESH_TOL)
    return INFTY;

  if (D >= 0.0) {
    D = std::sqrt(D);
    // the solution -p - D is always smaller as -p + D : Check this one first
    if (-p - D > l)
      return -p - D;
    if (-p + D > l)
      return -p + D;
  }

  return INFTY;
}

double SphericalMesh::find_theta_crossing(
  const Position& r, const Direction& u, double l, int shell) const
{
  // Theta grid is [0, π], thus there is no real surface to cross
  if (full_theta_ && (shape_[1] == 1))
    return INFTY;

  shell = sanitize_theta(shell);

  // solving z(s) = cos/theta) * r(s) with r(s) = r+s*u
  // yields
  // a*s^2 + 2*b*s + c == 0 with
  // a = cos(theta)^2 - u.z * u.z
  // b = r*u * cos(theta)^2 - u.z * r.z
  // c = r*r * cos(theta)^2 - r.z^2

  const double cos_t = std::cos(grid_[1][shell]);
  const bool sgn = std::signbit(cos_t);
  const double cos_t_2 = cos_t * cos_t;

  const double a = cos_t_2 - u.z * u.z;
  const double b = r.dot(u) * cos_t_2 - r.z * u.z;
  const double c = r.dot(r) * cos_t_2 - r.z * r.z;

  // if factor of s^2 is zero, direction of flight is parallel to theta surface
  if (std::abs(a) < FP_PRECISION) {
    // if b vanishes, direction of flight is within theta surface and crossing
    // is not possible
    if (std::abs(b) < FP_PRECISION)
      return INFTY;

    const double s = -0.5 * c / b;
    // Check if solution is in positive direction of flight and has correct sign
    if ((s > l) && (std::signbit(r.z + s * u.z) == sgn))
      return s;

    // no crossing is possible
    return INFTY;
  }

  const double p = b / a;
  double D = p * p - c / a;

  if (D < 0.0)
    return INFTY;

  D = std::sqrt(D);

  // the solution -p-D is always smaller as -p+D : Check this one first
  double s = -p - D;
  // Check if solution is in positive direction of flight and has correct sign
  if ((s > l) && (std::signbit(r.z + s * u.z) == sgn))
    return s;

  s = -p + D;
  // Check if solution is in positive direction of flight and has correct sign
  if ((s > l) && (std::signbit(r.z + s * u.z) == sgn))
    return s;

  return INFTY;
}

double SphericalMesh::find_phi_crossing(
  const Position& r, const Direction& u, double l, int shell) const
{
  // Phi grid is [0, 2π], thus there is no real surface to cross
  if (full_phi_ && (shape_[2] == 1))
    return INFTY;

  shell = sanitize_phi(shell);

  const double p0 = grid_[2][shell];

  // solve y(s)/x(s) = tan(p0) = sin(p0)/cos(p0)
  // => x(s) * cos(p0) = y(s) * sin(p0)
  // => (y + s * v) * cos(p0) = (x + s * u) * sin(p0)
  // = s * (v * cos(p0) - u * sin(p0)) = - (y * cos(p0) - x * sin(p0))

  const double c0 = std::cos(p0);
  const double s0 = std::sin(p0);

  const double denominator = (u.x * s0 - u.y * c0);

  // Check if direction of flight is not parallel to phi surface
  if (std::abs(denominator) > FP_PRECISION) {
    const double s = -(r.x * s0 - r.y * c0) / denominator;
    // Check if solution is in positive direction of flight and crosses the
    // correct phi surface (not -phi)
    if ((s > l) && ((c0 * (r.x + s * u.x) + s0 * (r.y + s * u.y)) > 0.0))
      return s;
  }

  return INFTY;
}

StructuredMesh::MeshDistance SphericalMesh::distance_to_grid_boundary(
  const MeshIndex& ijk, int i, const Position& r0, const Direction& u,
  double l) const
{

  if (i == 0) {
    return std::min(
      MeshDistance(ijk[i] + 1, true, find_r_crossing(r0, u, l, ijk[i])),
      MeshDistance(ijk[i] - 1, false, find_r_crossing(r0, u, l, ijk[i] - 1)));

  } else if (i == 1) {
    return std::min(MeshDistance(sanitize_theta(ijk[i] + 1), true,
                      find_theta_crossing(r0, u, l, ijk[i])),
      MeshDistance(sanitize_theta(ijk[i] - 1), false,
        find_theta_crossing(r0, u, l, ijk[i] - 1)));

  } else {
    return std::min(MeshDistance(sanitize_phi(ijk[i] + 1), true,
                      find_phi_crossing(r0, u, l, ijk[i])),
      MeshDistance(sanitize_phi(ijk[i] - 1), false,
        find_phi_crossing(r0, u, l, ijk[i] - 1)));
  }
}

int SphericalMesh::set_grid()
{
  shape_ = {static_cast<int>(grid_[0].size()) - 1,
    static_cast<int>(grid_[1].size()) - 1,
    static_cast<int>(grid_[2].size()) - 1};

  for (const auto& g : grid_) {
    if (g.size() < 2) {
      set_errmsg("x-, y-, and z- grids for spherical meshes "
                 "must each have at least 2 points");
      return OPENMC_E_INVALID_ARGUMENT;
    }
    if (std::adjacent_find(g.begin(), g.end(), std::greater_equal<>()) !=
        g.end()) {
      set_errmsg("Values in for r-, theta-, and phi- grids for "
                 "spherical meshes must be sorted and unique.");
      return OPENMC_E_INVALID_ARGUMENT;
    }
    if (g.front() < 0.0) {
      set_errmsg("r-, theta-, and phi- grids for "
                 "spherical meshes must start at v >= 0.");
      return OPENMC_E_INVALID_ARGUMENT;
    }
  }
  if (grid_[1].back() > PI) {
    set_errmsg("theta-grids for "
               "spherical meshes must end with theta <= pi.");

    return OPENMC_E_INVALID_ARGUMENT;
  }
  if (grid_[2].back() > 2 * PI) {
    set_errmsg("phi-grids for "
               "spherical meshes must end with phi <= 2*pi.");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  full_theta_ = (grid_[1].front() == 0.0) && (grid_[1].back() == PI);
  full_phi_ = (grid_[2].front() == 0.0) && (grid_[2].back() == 2 * PI);

  double r = grid_[0].back();
  lower_left_ = {origin_[0] - r, origin_[1] - r, origin_[2] - r};
  upper_right_ = {origin_[0] + r, origin_[1] + r, origin_[2] + r};

  return 0;
}

int SphericalMesh::get_index_in_direction(double r, int i) const
{
  return lower_bound_index(grid_[i].begin(), grid_[i].end(), r) + 1;
}

std::pair<vector<double>, vector<double>> SphericalMesh::plot(
  Position plot_ll, Position plot_ur) const
{
  fatal_error("Plot of spherical Mesh not implemented");

  // Figure out which axes lie in the plane of the plot.
  array<vector<double>, 2> axis_lines;
  return {axis_lines[0], axis_lines[1]};
}

void SphericalMesh::to_hdf5(hid_t group) const
{
  hid_t mesh_group = create_group(group, "mesh " + std::to_string(id_));

  write_dataset(mesh_group, "type", SphericalMesh::mesh_type);
  write_dataset(mesh_group, "r_grid", grid_[0]);
  write_dataset(mesh_group, "theta_grid", grid_[1]);
  write_dataset(mesh_group, "phi_grid", grid_[2]);
  write_dataset(mesh_group, "origin", origin_);

  close_group(mesh_group);
}

double SphericalMesh::volume(const MeshIndex& ijk) const
{
  double r_i = grid_[0][ijk[0] - 1];
  double r_o = grid_[0][ijk[0]];

  double theta_i = grid_[1][ijk[1] - 1];
  double theta_o = grid_[1][ijk[1]];

  double phi_i = grid_[2][ijk[2] - 1];
  double phi_o = grid_[2][ijk[2]];

  return (1.0 / 3.0) * (r_o * r_o * r_o - r_i * r_i * r_i) *
         (std::cos(theta_i) - std::cos(theta_o)) * (phi_o - phi_i);
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

template<class T>
bool is_mesh_type(int32_t index)
{
  T* mesh = dynamic_cast<T*>(model::meshes[index].get());
  return mesh;
}

//==============================================================================
// C API functions
//==============================================================================

// Return the type of mesh as a C string
extern "C" int openmc_mesh_get_type(int32_t index, char* type)
{
  if (int err = check_mesh(index))
    return err;

  std::strcpy(type, model::meshes[index].get()->get_mesh_type().c_str());

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
    if (RegularMesh::mesh_type == type) {
      model::meshes.push_back(make_unique<RegularMesh>());
    } else if (RectilinearMesh::mesh_type == type) {
      model::meshes.push_back(make_unique<RectilinearMesh>());
    } else if (CylindricalMesh::mesh_type == type) {
      model::meshes.push_back(make_unique<CylindricalMesh>());
    } else if (SphericalMesh::mesh_type == type) {
      model::meshes.push_back(make_unique<SphericalMesh>());
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
  if (lib_name == MOABMesh::mesh_lib_type) {
    model::meshes.push_back(std::move(make_unique<MOABMesh>(mesh_file)));
    valid_lib = true;
  }
#endif

#ifdef LIBMESH
  if (lib_name == LibMesh::mesh_lib_type) {
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

//! Get the number of elements in a mesh
extern "C" int openmc_mesh_get_n_elements(int32_t index, size_t* n)
{
  if (int err = check_mesh(index))
    return err;
  *n = model::meshes[index]->n_bins();
  return 0;
}

//! Get the volume of each element in the mesh
extern "C" int openmc_mesh_get_volumes(int32_t index, double* volumes)
{
  if (int err = check_mesh(index))
    return err;
  for (int i = 0; i < model::meshes[index]->n_bins(); ++i) {
    volumes[i] = model::meshes[index]->volume(i);
  }
  return 0;
}

//! Get the bounding box of a mesh
extern "C" int openmc_mesh_bounding_box(int32_t index, double* ll, double* ur)
{
  if (int err = check_mesh(index))
    return err;

  BoundingBox bbox = model::meshes[index]->bounding_box();

  // set lower left corner values
  ll[0] = bbox.xmin;
  ll[1] = bbox.ymin;
  ll[2] = bbox.zmin;

  // set upper right corner values
  ur[0] = bbox.xmax;
  ur[1] = bbox.ymax;
  ur[2] = bbox.zmax;
  return 0;
}

extern "C" int openmc_mesh_material_volumes(int32_t index, int n_sample,
  int bin, int result_size, void* result, int* hits, uint64_t* seed)
{
  auto result_ = reinterpret_cast<Mesh::MaterialVolume*>(result);
  if (!result_) {
    set_errmsg("Invalid result pointer passed to openmc_mesh_material_volumes");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  if (int err = check_mesh(index))
    return err;

  int n = model::meshes[index]->material_volumes(
    n_sample, bin, {result_, result_ + result_size}, seed);
  *hits = n;
  return (n == -1) ? OPENMC_E_ALLOCATE : 0;
}

extern "C" int openmc_mesh_get_plot_bins(int32_t index, Position origin,
  Position width, int basis, int* pixels, int32_t* data)
{
  if (int err = check_mesh(index))
    return err;
  const auto& mesh = model::meshes[index].get();

  int pixel_width = pixels[0];
  int pixel_height = pixels[1];

  // get pixel size
  double in_pixel = (width[0]) / static_cast<double>(pixel_width);
  double out_pixel = (width[1]) / static_cast<double>(pixel_height);

  // setup basis indices and initial position centered on pixel
  int in_i, out_i;
  Position xyz = origin;
  enum class PlotBasis { xy = 1, xz = 2, yz = 3 };
  PlotBasis basis_enum = static_cast<PlotBasis>(basis);
  switch (basis_enum) {
  case PlotBasis::xy:
    in_i = 0;
    out_i = 1;
    break;
  case PlotBasis::xz:
    in_i = 0;
    out_i = 2;
    break;
  case PlotBasis::yz:
    in_i = 1;
    out_i = 2;
    break;
  default:
    UNREACHABLE();
  }

  // set initial position
  xyz[in_i] = origin[in_i] - width[0] / 2. + in_pixel / 2.;
  xyz[out_i] = origin[out_i] + width[1] / 2. - out_pixel / 2.;

#pragma omp parallel
  {
    Position r = xyz;

#pragma omp for
    for (int y = 0; y < pixel_height; y++) {
      r[out_i] = xyz[out_i] - out_pixel * y;
      for (int x = 0; x < pixel_width; x++) {
        r[in_i] = xyz[in_i] + in_pixel * x;
        data[pixel_width * y + x] = mesh->get_bin(r);
      }
    }
  }

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
  mesh->n_dimension_ = n;
  std::copy(dims, dims + n, mesh->shape_.begin());
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

  // Set material volumes

  // TODO: incorporate this into method in RegularMesh that can be called from
  // here and from constructor
  m->volume_frac_ = 1.0 / xt::prod(m->get_x_shape())();
  m->element_volume_ = 1.0;
  for (int i = 0; i < m->n_dimension_; i++) {
    m->element_volume_ *= m->width_[i];
  }

  return 0;
}

//! Set the mesh parameters for rectilinear, cylindrical and spharical meshes
template<class C>
int openmc_structured_mesh_set_grid_impl(int32_t index, const double* grid_x,
  const int nx, const double* grid_y, const int ny, const double* grid_z,
  const int nz)
{
  if (int err = check_mesh_type<C>(index))
    return err;

  C* m = dynamic_cast<C*>(model::meshes[index].get());

  m->n_dimension_ = 3;

  m->grid_[0].reserve(nx);
  m->grid_[1].reserve(ny);
  m->grid_[2].reserve(nz);

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

//! Get the mesh parameters for rectilinear, cylindrical and spherical meshes
template<class C>
int openmc_structured_mesh_get_grid_impl(int32_t index, double** grid_x,
  int* nx, double** grid_y, int* ny, double** grid_z, int* nz)
{
  if (int err = check_mesh_type<C>(index))
    return err;
  C* m = dynamic_cast<C*>(model::meshes[index].get());

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

//! Get the rectilinear mesh grid
extern "C" int openmc_rectilinear_mesh_get_grid(int32_t index, double** grid_x,
  int* nx, double** grid_y, int* ny, double** grid_z, int* nz)
{
  return openmc_structured_mesh_get_grid_impl<RectilinearMesh>(
    index, grid_x, nx, grid_y, ny, grid_z, nz);
}

//! Set the rectilienar mesh parameters
extern "C" int openmc_rectilinear_mesh_set_grid(int32_t index,
  const double* grid_x, const int nx, const double* grid_y, const int ny,
  const double* grid_z, const int nz)
{
  return openmc_structured_mesh_set_grid_impl<RectilinearMesh>(
    index, grid_x, nx, grid_y, ny, grid_z, nz);
}

//! Get the cylindrical mesh grid
extern "C" int openmc_cylindrical_mesh_get_grid(int32_t index, double** grid_x,
  int* nx, double** grid_y, int* ny, double** grid_z, int* nz)
{
  return openmc_structured_mesh_get_grid_impl<CylindricalMesh>(
    index, grid_x, nx, grid_y, ny, grid_z, nz);
}

//! Set the cylindrical mesh parameters
extern "C" int openmc_cylindrical_mesh_set_grid(int32_t index,
  const double* grid_x, const int nx, const double* grid_y, const int ny,
  const double* grid_z, const int nz)
{
  return openmc_structured_mesh_set_grid_impl<CylindricalMesh>(
    index, grid_x, nx, grid_y, ny, grid_z, nz);
}

//! Get the spherical mesh grid
extern "C" int openmc_spherical_mesh_get_grid(int32_t index, double** grid_x,
  int* nx, double** grid_y, int* ny, double** grid_z, int* nz)
{

  return openmc_structured_mesh_get_grid_impl<SphericalMesh>(
    index, grid_x, nx, grid_y, ny, grid_z, nz);
  ;
}

//! Set the spherical mesh parameters
extern "C" int openmc_spherical_mesh_set_grid(int32_t index,
  const double* grid_x, const int nx, const double* grid_y, const int ny,
  const double* grid_z, const int nz)
{
  return openmc_structured_mesh_set_grid_impl<SphericalMesh>(
    index, grid_x, nx, grid_y, ny, grid_z, nz);
}

#ifdef DAGMC

const std::string MOABMesh::mesh_lib_type = "moab";

MOABMesh::MOABMesh(pugi::xml_node node) : UnstructuredMesh(node)
{
  initialize();
}

MOABMesh::MOABMesh(const std::string& filename, double length_multiplier)
{
  filename_ = filename;
  set_length_multiplier(length_multiplier);
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

  // set member range of vertices
  int vertex_dim = 0;
  rval = mbi_->get_entities_by_dimension(0, vertex_dim, verts_);
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to get all vertex handles");
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

  if (length_multiplier_ > 0.0) {
    // get the connectivity of all tets
    moab::Range adj;
    rval = mbi_->get_adjacencies(ehs_, 0, true, adj, moab::Interface::UNION);
    if (rval != moab::MB_SUCCESS) {
      fatal_error("Failed to get adjacent vertices of tetrahedra.");
    }
    // scale all vertex coords by multiplier (done individually so not all
    // coordinates are in memory twice at once)
    for (auto vert : adj) {
      // retrieve coords
      std::array<double, 3> coord;
      rval = mbi_->get_coords(&vert, 1, coord.data());
      if (rval != moab::MB_SUCCESS) {
        fatal_error("Could not get coordinates of vertex.");
      }
      // scale coords
      for (auto& c : coord) {
        c *= length_multiplier_;
      }
      // set new coords
      rval = mbi_->set_coords(&vert, 1, coord.data());
      if (rval != moab::MB_SUCCESS) {
        fatal_error("Failed to set new vertex coordinates");
      }
    }
  }

  // Determine bounds of mesh
  this->determine_bounds();
}

void MOABMesh::prepare_for_tallies()
{
  // if the KDTree has already been constructed, do nothing
  if (kdtree_)
    return;

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
  write_message("Getting tet adjacencies...", 7);
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
  write_message("Building adaptive k-d tree for tet mesh...", 7);
  kdtree_ = make_unique<moab::AdaptiveKDTree>(mbi_.get());

  // Determine what options to use
  std::ostringstream options_stream;
  if (options_.empty()) {
    options_stream << "MAX_DEPTH=20;PLANE_SET=2;";
  } else {
    options_stream << options_;
  }
  moab::FileOptions file_opts(options_stream.str().c_str());

  // Build the k-d tree
  rval = kdtree_->build_tree(all_tets_and_tris, &kdtree_root_, &file_opts);
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
  return mesh_lib_type;
}

// Sample position within a tet for MOAB type tets
Position MOABMesh::sample_element(int32_t bin, uint64_t* seed) const
{

  moab::EntityHandle tet_ent = get_ent_handle_from_bin(bin);

  // Get vertex coordinates for MOAB tet
  const moab::EntityHandle* conn1;
  int conn1_size;
  moab::ErrorCode rval = mbi_->get_connectivity(tet_ent, conn1, conn1_size);
  if (rval != moab::MB_SUCCESS || conn1_size != 4) {
    fatal_error(fmt::format(
      "Failed to get tet connectivity or connectivity size ({}) is invalid.",
      conn1_size));
  }
  moab::CartVect p[4];
  rval = mbi_->get_coords(conn1, conn1_size, p[0].array());
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to get tet coords");
  }

  std::array<Position, 4> tet_verts;
  for (int i = 0; i < 4; i++) {
    tet_verts[i] = {p[i][0], p[i][1], p[i][2]};
  }
  // Samples position within tet using Barycentric stuff
  return this->sample_tet(tet_verts, seed);
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

int MOABMesh::get_vert_idx_from_handle(moab::EntityHandle vert) const
{
  int idx = vert - verts_[0];
  if (idx >= n_vertices()) {
    fatal_error(
      fmt::format("Invalid vertex idx {} (# vertices {})", idx, n_vertices()));
  }
  return idx;
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

int MOABMesh::n_vertices() const
{
  return verts_.size();
}

Position MOABMesh::vertex(int id) const
{

  moab::ErrorCode rval;

  moab::EntityHandle vert = verts_[id];

  moab::CartVect coords;
  rval = mbi_->get_coords(&vert, 1, coords.array());
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to get the coordinates of a vertex.");
  }

  return {coords[0], coords[1], coords[2]};
}

std::vector<int> MOABMesh::connectivity(int bin) const
{
  moab::ErrorCode rval;

  auto tet = get_ent_handle_from_bin(bin);

  // look up the tet connectivity
  vector<moab::EntityHandle> conn;
  rval = mbi_->get_connectivity(&tet, 1, conn);
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to get connectivity of a mesh element.");
    return {};
  }

  std::vector<int> verts(4);
  for (int i = 0; i < verts.size(); i++) {
    verts[i] = get_vert_idx_from_handle(conn[i]);
  }

  return verts;
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

const std::string LibMesh::mesh_lib_type = "libmesh";

LibMesh::LibMesh(pugi::xml_node node) : UnstructuredMesh(node)
{
  // filename_ and length_multiplier_ will already be set by the
  // UnstructuredMesh constructor
  set_mesh_pointer_from_filename(filename_);
  set_length_multiplier(length_multiplier_);
  initialize();
}

// create the mesh from a pointer to a libMesh Mesh
LibMesh::LibMesh(libMesh::MeshBase& input_mesh, double length_multiplier)
{
  m_ = &input_mesh;
  set_length_multiplier(length_multiplier);
  initialize();
}

// create the mesh from an input file
LibMesh::LibMesh(const std::string& filename, double length_multiplier)
{
  set_mesh_pointer_from_filename(filename);
  set_length_multiplier(length_multiplier);
  initialize();
}

void LibMesh::set_mesh_pointer_from_filename(const std::string& filename)
{
  filename_ = filename;
  unique_m_ = make_unique<libMesh::Mesh>(*settings::libmesh_comm, n_dimension_);
  m_ = unique_m_.get();
  m_->read(filename_);
}

// intialize from mesh file
void LibMesh::initialize()
{
  if (!settings::libmesh_comm) {
    fatal_error(
      "Attempting to use an unstructured mesh without a libMesh communicator.");
  }

  // assuming that unstructured meshes used in OpenMC are 3D
  n_dimension_ = 3;

  if (length_multiplier_ > 0.0) {
    libMesh::MeshTools::Modification::scale(*m_, length_multiplier_);
  }
  // if OpenMC is managing the libMesh::MeshBase instance, prepare the mesh.
  // Otherwise assume that it is prepared by its owning application
  if (unique_m_) {
    m_->prepare_for_use();
  }

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

  for (int i = 0; i < num_threads(); i++) {
    pl_.emplace_back(m_->sub_point_locator());
    pl_.back()->set_contains_point_tol(FP_COINCIDENT);
    pl_.back()->enable_out_of_mesh_mode();
  }

  // store first element in the mesh to use as an offset for bin indices
  auto first_elem = *m_->elements_begin();
  first_element_id_ = first_elem->id();

  // bounding box for the mesh for quick rejection checks
  bbox_ = libMesh::MeshTools::create_bounding_box(*m_);
  libMesh::Point ll = bbox_.min();
  libMesh::Point ur = bbox_.max();
  lower_left_ = {ll(0), ll(1), ll(2)};
  upper_right_ = {ur(0), ur(1), ur(2)};
}

// Sample position within a tet for LibMesh type tets
Position LibMesh::sample_element(int32_t bin, uint64_t* seed) const
{
  const auto& elem = get_element_from_bin(bin);
  // Get tet vertex coordinates from LibMesh
  std::array<Position, 4> tet_verts;
  for (int i = 0; i < elem.n_nodes(); i++) {
    auto node_ref = elem.node_ref(i);
    tet_verts[i] = {node_ref(0), node_ref(1), node_ref(2)};
  }
  // Samples position within tet using Barycentric coordinates
  return this->sample_tet(tet_verts, seed);
}

Position LibMesh::centroid(int bin) const
{
  const auto& elem = this->get_element_from_bin(bin);
  auto centroid = elem.vertex_average();
  return {centroid(0), centroid(1), centroid(2)};
}

int LibMesh::n_vertices() const
{
  return m_->n_nodes();
}

Position LibMesh::vertex(int vertex_id) const
{
  const auto node_ref = m_->node_ref(vertex_id);
  return {node_ref(0), node_ref(1), node_ref(2)};
}

std::vector<int> LibMesh::connectivity(int elem_id) const
{
  std::vector<int> conn;
  const auto* elem_ptr = m_->elem_ptr(elem_id);
  for (int i = 0; i < elem_ptr->n_nodes(); i++) {
    conn.push_back(elem_ptr->node_id(i));
  }
  return conn;
}

std::string LibMesh::library() const
{
  return mesh_lib_type;
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

  const auto& point_locator = pl_.at(thread_num());

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
  return this->get_element_from_bin(bin).volume();
}

#endif // LIBMESH

//==============================================================================
// Non-member functions
//==============================================================================

void read_meshes(pugi::xml_node root)
{
  std::unordered_set<int> mesh_ids;

  for (auto node : root.children("mesh")) {
    // Check to make sure multiple meshes in the same file don't share IDs
    int id = std::stoi(get_node_value(node, "id"));
    if (contains(mesh_ids, id)) {
      fatal_error(fmt::format(
        "Two or more meshes use the same unique ID '{}' in the same input file",
        id));
    }
    mesh_ids.insert(id);

    // If we've already read a mesh with the same ID in a *different* file,
    // assume it is the same here
    if (model::mesh_map.find(id) != model::mesh_map.end()) {
      warning(fmt::format("Mesh with ID={} appears in multiple files.", id));
      continue;
    }

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
    if (mesh_type == RegularMesh::mesh_type) {
      model::meshes.push_back(make_unique<RegularMesh>(node));
    } else if (mesh_type == RectilinearMesh::mesh_type) {
      model::meshes.push_back(make_unique<RectilinearMesh>(node));
    } else if (mesh_type == CylindricalMesh::mesh_type) {
      model::meshes.push_back(make_unique<CylindricalMesh>(node));
    } else if (mesh_type == SphericalMesh::mesh_type) {
      model::meshes.push_back(make_unique<SphericalMesh>(node));
#ifdef DAGMC
    } else if (mesh_type == UnstructuredMesh::mesh_type &&
               mesh_lib == MOABMesh::mesh_lib_type) {
      model::meshes.push_back(make_unique<MOABMesh>(node));
#endif
#ifdef LIBMESH
    } else if (mesh_type == UnstructuredMesh::mesh_type &&
               mesh_lib == LibMesh::mesh_lib_type) {
      model::meshes.push_back(make_unique<LibMesh>(node));
#endif
    } else if (mesh_type == UnstructuredMesh::mesh_type) {
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
