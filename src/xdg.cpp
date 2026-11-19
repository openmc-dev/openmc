#include "openmc/xdg.h"

#include <algorithm> // for remove_if
#include <cmath>     // for abs
#include <memory>    // for shared_ptr
#include <string>    // for string
#include <utility>   // for pair
#include <vector>    // for vector

#include <fmt/core.h>

#ifdef OPENMC_XDG_ENABLED
#include "xdg/config.h"
#include "xdg/xdg.h"
#endif

#include "openmc/constants.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/geometry_aux.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/settings.h"
#include "openmc/string_utils.h"

namespace openmc {

#ifdef OPENMC_XDG_ENABLED
const bool XDG_ENABLED = true;
#else
const bool XDG_ENABLED = false;
#endif

} // namespace openmc

#ifdef OPENMC_XDG_ENABLED

namespace openmc {

//==============================================================================
// XDG Mesh implementation
//==============================================================================

const std::string XDGMesh::mesh_interface = "xdg";

XDGMesh::XDGMesh(pugi::xml_node node) : UnstructuredMesh(node)
{
  std::string mesh_lib = get_node_value(node, "library", true, true);
  if (mesh_lib == "moab") {
    mesh_library_ = xdg::MeshLibrary::MOAB;
  } else if (mesh_lib == "libmesh") {
    mesh_library_ = xdg::MeshLibrary::LIBMESH;
  } else {
    fatal_error(fmt::format("Invalid mesh library specified in XML: {}. "
                            "Valid options are 'moab' and 'libmesh'.",
      mesh_lib));
  }

  initialize();
}

XDGMesh::XDGMesh(hid_t group) : UnstructuredMesh(group)
{
  std::string mesh_lib;
  read_dataset(group, "library", mesh_lib);
  if (mesh_lib == "moab") {
    mesh_library_ = xdg::MeshLibrary::MOAB;
  } else if (mesh_lib == "libmesh") {
    mesh_library_ = xdg::MeshLibrary::LIBMESH;
  } else {
    fatal_error(fmt::format("Invalid mesh library specified in HDF5: {}. "
                            "Valid options are 'moab' and 'libmesh'.",
      mesh_lib));
  }
  initialize();
}

XDGMesh::XDGMesh(const std::string& filename, double length_multiplier)
{
  filename_ = filename;
  set_length_multiplier(length_multiplier);
  initialize();
}

XDGMesh::XDGMesh(
  std::shared_ptr<xdg::XDG> external_xdg, double length_multiplier)
{
  xdg_ = external_xdg;
  filename_ = "unknown (external file)";
  set_length_multiplier(length_multiplier);
  initialize();
}

void XDGMesh::initialize()
{
  interface_ = "xdg";

  // Externally supplied instances already have an initialized mesh manager.
  if (!xdg_) {
    if (!xdg::XDGConfig::config().mesh_manager_enabled(mesh_library_)) {
      fatal_error(fmt::format("Mesh library {} is not enabled in XDG.",
        mesh_library_ == xdg::MeshLibrary::LIBMESH ? "libMesh" : "MOAB"));
    }

    xdg_ = xdg::XDG::create(mesh_library_);
    if (!file_exists(filename_)) {
      fatal_error(fmt::format(
        "Mesh file \"{}\" for mesh {} does not exist", filename_, id_));
    }

    xdg_->mesh_manager()->load_file(filename_);
    xdg_->mesh_manager()->init();
    xdg_->mesh_manager()->parse_metadata();

    // A file with only surfaces, such as a DAGMC geometry, loads without error
    if (xdg_->mesh_manager()->num_volume_elements() == 0) {
      fatal_error(fmt::format(
        "Mesh file \"{}\" for mesh {} has no volume elements", filename_, id_));
    }
  }

  // Cache the number of volume element surfaces.
  n_surface_bins_ = xdg_->mesh_manager()->num_volume_element_faces();

  auto global_bbox = xdg_->mesh_manager()->global_bounding_box();
  double length_multiplier =
    length_multiplier_ > 0.0 ? length_multiplier_ : 1.0;
  Position lower_left =
    length_multiplier *
    Position(global_bbox.min_x, global_bbox.min_y, global_bbox.min_z);
  Position upper_right =
    length_multiplier *
    Position(global_bbox.max_x, global_bbox.max_y, global_bbox.max_z);
  lower_left_ = {lower_left.x, lower_left.y, lower_left.z};
  upper_right_ = {upper_right.x, upper_right.y, upper_right.z};
}

void XDGMesh::prepare_for_point_location()
{
  // A mesh can be used by several mesh filters, weight windows and random ray
  // subdivision, which each prepare it; preparing XDG again would register
  // its volumes and build its trees again
  if (prepared_)
    return;
  xdg_->prepare_raytracer();
  prepared_ = true;
}

Position XDGMesh::sample_element(int32_t bin, uint64_t* seed) const
{
  auto vertices = xdg_->mesh_manager()->element_vertices(bin_to_mesh_id(bin));
  if (vertices.size() != 4) {
    fatal_error("Sampling from an XDG mesh is only supported for linear "
                "tetrahedra.");
  }
  Position sampled_position = this->sample_tet<xdg::Vertex>(vertices, seed);
  if (length_multiplier_ > 0.0) {
    sampled_position *= length_multiplier_;
  }
  return sampled_position;
}

void XDGMesh::bins_crossed(Position r0, Position r1, const Direction& u,
  vector<int>& bins, vector<double>& lengths) const
{
  xdg::Position p0 {r0.x, r0.y, r0.z};
  xdg::Position p1 {r1.x, r1.y, r1.z};
  if (length_multiplier_ > 0.0) {
    p0 /= length_multiplier_;
    p1 /= length_multiplier_;
  }
  double inv_length = 1 / (p1 - p0).length();
  auto track_segments = xdg_->segments(p0, p1);
  // remove elements with lengths of zero
  track_segments.erase(
    std::remove_if(track_segments.begin(), track_segments.end(),
      [](const std::pair<xdg::MeshID, double>& p) { return p.second == 0.0; }),
    track_segments.end());
  for (const auto& track_segment : track_segments) {
    bins.push_back(mesh_id_to_bin(track_segment.first));
    lengths.push_back(track_segment.second * inv_length);
  }
}

//! Split the track from r0 to r1 into the mesh elements it crosses and the
//! stretches outside of the mesh (bin C_NONE), in order along the track, with
//! lengths as fractions of the track's length. The track is walked from its
//! start, repeating two steps until its end is reached:
//!  - If the current point is in no element, fire a ray against the mesh's
//!    implicit complement to find where the track next enters the mesh. The
//!    stretch up to there is outside of the mesh, as is the rest of the track
//!    if it never enters.
//!  - From the element there, walk the elements along the track until it
//!    leaves the mesh or ends, recording each element crossed.
//! The element that the last walk ended in is not taken again at the next
//! lookup, so a point on its exit face, or a walk that got nowhere, does not
//! repeat forever.
void XDGMesh::bins_crossed_with_outside(Position r0, Position r1,
  const Direction& u, vector<int>& bins, vector<double>& lengths) const
{
  xdg::Position r {r0.x, r0.y, r0.z};
  xdg::Position r_end {r1.x, r1.y, r1.z};
  if (length_multiplier_ > 0.0) {
    r /= length_multiplier_;
    r_end /= length_multiplier_;
  }
  xdg::Direction direction = r_end - r;
  const double total = direction.length();
  if (total == 0.0)
    return;
  direction /= total;

  const auto add = [&](int bin, double length) {
    bins.push_back(bin);
    lengths.push_back(length / total);
  };

  // This is the loop of xdg::XDG::segments, which walks the elements from where
  // the track enters the mesh each time, but does not say how far it skipped to
  // get there. It can go once XDG does.
  const auto& manager = xdg_->mesh_manager();
  const xdg::MeshID complement = manager->implicit_complement();
  std::vector<xdg::MeshID> hit_faces;
  double distance = total;
  xdg::MeshID last = xdg::ID_NONE;
  do {
    // A lookup on the face that the last walk left through can give that
    // element again, as can one where a walk did not get anywhere
    xdg::MeshID element = xdg_->find_element(r);
    if (element == last)
      element = xdg::ID_NONE;

    if (element == xdg::ID_NONE) {
      const auto hit = xdg_->ray_fire(complement, r + direction * xdg::TINY_BIT,
        direction, xdg::INFTY, xdg::HitOrientation::EXITING, &hit_faces);
      const double to_mesh = hit.first + xdg::TINY_BIT;
      if (hit.second == xdg::ID_NONE || to_mesh >= distance) {
        add(C_NONE, distance);
        return;
      }
      add(C_NONE, to_mesh);
      r += direction * to_mesh;
      distance -= to_mesh;
      element = manager->get_boundary_face_element(hit_faces.back());
      hit_faces.clear();
      if (element == xdg::ID_NONE) {
        warning(fmt::format(
          "A track enters mesh {} through a face with no element.", id_));
        add(C_NONE, distance);
        return;
      }
    }

    const auto walk = manager->walk_elements(element, r, direction, distance);
    double walked = 0.0;
    for (const auto& [id, length] : walk) {
      if (length > 0.0) {
        add(mesh_id_to_bin(id), length);
        walked += length;
      }
    }
    last = walked > 0.0 ? walk.back().first : element;
    r += direction * walked;
    distance -= walked;
    // (After a walk to the end of the track, what is left is roundoff)
  } while (distance > xdg::TINY_BIT);
}

int XDGMesh::get_bin(Position r) const
{
  xdg::Position p {r.x, r.y, r.z};
  if (length_multiplier_ > 0.0) {
    p /= length_multiplier_;
  }
  return mesh_id_to_bin(xdg_->find_element(p));
}

int XDGMesh::n_bins() const
{
  return xdg_->mesh_manager()->num_volume_elements();
}

int XDGMesh::n_surface_bins() const
{
  return n_surface_bins_;
}

std::pair<vector<double>, vector<double>> XDGMesh::plot(
  Position plot_ll, Position plot_ur) const
{
  fatal_error("Plot of XDGMesh mesh not implemented");
  return {};
}

std::string XDGMesh::library() const
{
  if (mesh_library_ == xdg::MeshLibrary::LIBMESH) {
    return "libmesh";
  } else if (mesh_library_ == xdg::MeshLibrary::MOAB) {
    return "moab";
  }

  fatal_error("Invalid mesh library specified for XDGMesh.");
  return {};
}

void XDGMesh::write(const std::string& base_filename) const
{
  warning("XDGMesh mesh write from C++ not implemented");
}

Position XDGMesh::centroid(int bin) const
{
  auto element_vertices =
    xdg_->mesh_manager()->element_vertices(bin_to_mesh_id(bin));

  xdg::Vertex center {0.0, 0.0, 0.0};
  for (const auto& v : element_vertices) {
    center += v;
  }

  center /= double(element_vertices.size());

  if (length_multiplier_ > 0.0) {
    center *= length_multiplier_;
  }

  return {center[0], center[1], center[2]};
}

int XDGMesh::n_vertices() const
{
  return xdg_->mesh_manager()->num_vertices();
}

Position XDGMesh::vertex(int id) const
{
  xdg::MeshID mesh_id = xdg_->mesh_manager()->vertex_id(id);
  auto v = xdg_->mesh_manager()->vertex_coordinates(mesh_id);
  if (length_multiplier_ > 0.0) {
    v *= length_multiplier_;
  }
  return {v[0], v[1], v[2]};
}

std::vector<int> XDGMesh::connectivity(int id) const
{
  auto conn = xdg_->mesh_manager()->element_connectivity(bin_to_mesh_id(id));
  for (auto& c : conn) {
    c = xdg_->mesh_manager()->vertex_index(c);
  }
  return conn;
}

double XDGMesh::volume(int bin) const
{
  // Positive whatever the element's node order
  double vol =
    std::abs(xdg_->mesh_manager()->element_volume(bin_to_mesh_id(bin)));
  if (length_multiplier_ > 0.0) {
    vol *= length_multiplier_ * length_multiplier_ * length_multiplier_;
  }
  return vol;
}

xdg::MeshID XDGMesh::bin_to_mesh_id(int bin) const
{
  return xdg_->mesh_manager()->element_id(bin);
}

int32_t XDGMesh::mesh_id_to_bin(xdg::MeshID id) const
{
  if (id < 0)
    return -1;
  return xdg_->mesh_manager()->element_index(id);
}

} // namespace openmc

#endif // XDG
