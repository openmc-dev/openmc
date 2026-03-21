#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

#include "openmc/xdg.h"

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

#include <fmt/core.h>

#ifdef OPENMC_XDG_ENABLED
#include "xdg/xdg.h"
#endif

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

const std::string XDGMesh::mesh_lib_type = "xdg";

XDGMesh::XDGMesh(pugi::xml_node node) : UnstructuredMesh(node)
{
  std::string mesh_lib = get_node_value(node, "library", true, true);
  if (mesh_lib == "moab") {
    mesh_library_ = xdg::MeshLibrary::MOAB;
  } else if (mesh_lib == "libmesh") {
    mesh_library_ = xdg::MeshLibrary::LIBMESH;
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
  }
  initialize();
}

XDGMesh::XDGMesh(const std::string& filename, double length_multiplier)
{
  filename_ = filename;
  set_length_multiplier(length_multiplier);
  initialize();
}

XDGMesh::XDGMesh(std::shared_ptr<xdg::XDG> external_xdg)
{
  xdg_ = external_xdg;
  filename_ = "unknown (external file)";
  initialize();
}

void XDGMesh::initialize()
{
  if (xdg_)
    return;

  // create XDGMesh instance
  xdg_ = xdg::XDG::create(mesh_library_);

  // load XDGMesh file
  if (!file_exists(filename_)) {
    fatal_error(fmt::format("Mesh file \"{}\" does not exist", filename_));
  }

  xdg_->mesh_manager()->load_file(filename_);
  xdg_->mesh_manager()->init();
  xdg_->mesh_manager()->parse_metadata();
}

void XDGMesh::prepare_for_point_location()
{
  xdg_->prepare_raytracer();
}

Position XDGMesh::sample_element(int32_t bin, uint64_t* seed) const
{
  // MeshIDs are 1-indexed, so we add 1 to the bin, which is 0-indexed
  auto vertices = xdg_->mesh_manager()->element_vertices(bin_to_mesh_id(bin));
  return this->sample_tet<xdg::Vertex>(vertices, seed);
}

void XDGMesh::bins_crossed(Position r0, Position r1, const Direction& u,
  vector<int>& bins, vector<double>& lengths) const
{
  // TODO: Make more robust (including mesh entrance/re-entrance)
  xdg::Position p0 {r0.x, r0.y, r0.z};
  xdg::Position p1 {r1.x, r1.y, r1.z};
  double length_rcp = 1 / (p1 - p0).length();
  auto track_segments = xdg_->segments(p0, p1);
  // remove elements with lengths of zero
  track_segments.erase(
    std::remove_if(track_segments.begin(), track_segments.end(),
      [](const std::pair<xdg::MeshID, double>& p) { return p.second == 0.0; }),
    track_segments.end());
  for (const auto& track_segment : track_segments) {
    bins.push_back(mesh_id_to_bin(track_segment.first));
    lengths.push_back(track_segment.second * length_rcp);
  }
}

int XDGMesh::get_bin(Position r) const
{
  xdg::Position p {r.x, r.y, r.z};
  return mesh_id_to_bin(xdg_->find_element(p));
}

int XDGMesh::n_bins() const
{
  return xdg_->mesh_manager()->num_volume_elements();
}

int XDGMesh::n_surface_bins() const
{
  return 4 * n_bins();
}

std::pair<vector<double>, vector<double>> XDGMesh::plot(
  Position plot_ll, Position plot_ur) const
{
  fatal_error("Plot of XDGMesh mesh not implemented");
  return {};
}

std::string XDGMesh::library() const
{
  return mesh_lib_type;
}

std::string XDGMesh::mesh_library() const
{
  if (mesh_library_ == xdg::MeshLibrary::LIBMESH) {
    return "libmesh";
  } else if (mesh_library_ == xdg::MeshLibrary::MOAB) {
    return "moab";
  }
}

void XDGMesh::write(const std::string& base_filename) const
{
  warning("XDGMesh mesh write from C++ not implemented");
}

Position XDGMesh::centroid(int bin) const
{
  auto element_vertices =
    xdg_->mesh_manager()->element_vertices(bin_to_mesh_id(bin));

  xdg::Vertex centroid {0.0, 0.0, 0.0};
  for (const auto& v : element_vertices) {
    centroid += v;
  }

  centroid /= double(element_vertices.size());

  return {centroid[0], centroid[1], centroid[2]};
}

int XDGMesh::n_vertices() const
{
  return xdg_->mesh_manager()->num_vertices();
}

Position XDGMesh::vertex(int id) const
{
  xdg::MeshID mesh_id = xdg_->mesh_manager()->vertex_id(id);
  auto v = xdg_->mesh_manager()->vertex_coordinates(mesh_id);
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
  return xdg_->mesh_manager()->element_volume(bin_to_mesh_id(bin));
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
