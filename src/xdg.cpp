#include "openmc/xdg.h"

#include <algorithm> // for remove_if
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
// XDG Universe implementation
//==============================================================================

XDGUniverse::XDGUniverse(pugi::xml_node node)
{
  if (check_for_node(node, "id")) {
    id_ = std::stoi(get_node_value(node, "id"));
  } else {
    fatal_error("Must specify the id of the XDG universe");
  }

  // TODO: fix this to use the mesh itself
  if (check_for_node(node, "mesh")) {
    int32_t mesh_id = std::stoi(get_node_value(node, "mesh"));
    const auto& mesh = model::meshes[model::mesh_map.at(mesh_id)];
    const auto& xdg_mesh = dynamic_cast<XDGMesh*>(mesh.get());
    if (xdg_mesh == nullptr) {
      fatal_error("XDG universe must be associated with an XDG mesh");
    }
    mesh_idx_ = model::mesh_map.at(mesh_id);
  }

  adjust_geometry_ids_ = false;
  if (check_for_node(node, "auto_geom_ids")) {
    adjust_geometry_ids_ = get_node_value_bool(node, "auto_geom_ids");
  }

  adjust_material_ids_ = false;
  if (check_for_node(node, "auto_mat_ids")) {
    adjust_material_ids_ = get_node_value_bool(node, "auto_mat_ids");
  }

  initialize();
}

// XDGUniverse::XDGUniverse(
//   const std::string& filename, bool auto_geom_ids, bool auto_mat_ids)
//   : filename_(filename), adjust_geometry_ids_(auto_geom_ids),
//     adjust_material_ids_(auto_mat_ids)
// {
//   geom_type() = GeometryType::XDG_SURFACE_MESH;
//   set_id();
//   initialize();
// }

// XDGUniverse::XDGUniverse(std::shared_ptr<xdg::XDG> xdg_ptr,
//   const std::string& filename, bool auto_geom_ids, bool auto_mat_ids)
//   : xdg_instance_(xdg_ptr), filename_(filename),
//     adjust_geometry_ids_(auto_geom_ids), adjust_material_ids_(auto_mat_ids)
// {
//   geom_type() = GeometryType::XDG_SURFACE_MESH;
//   set_id();
//   init_metadata();
//   init_geometry();
// }

void XDGUniverse::set_id()
{
  // determine the next universe id
  int32_t next_univ_id = 0;
  for (const auto& u : model::universes) {
    if (u->id_ > next_univ_id)
      next_univ_id = u->id_;
  }
  next_univ_id++;

  // set the universe id
  id_ = next_univ_id;
}

void XDGUniverse::initialize()
{
  init_geometry();
  mesh()->prepare_for_point_location();
}

void XDGUniverse::init_xdg()
{
  // create a new XDG instance
  if (library() == "moab")
    xdg_instance_ = xdg::XDG::create(xdg::MeshLibrary::MOAB);
  else if (library() == "libmesh")
    xdg_instance_ = xdg::XDG::create(xdg::MeshLibrary::LIBMESH);
  else
    fatal_error("Unknown XDG library specified: " + library());

    // load the XDG geometry
  if (!file_exists(filename())) {
    fatal_error("Geometry XDG file '" + filename() + "' does not exist!");
  }
  xdg_instance_->mesh_manager()->load_file(filename());
  xdg_instance_->mesh_manager()->init();
  xdg_instance_->prepare_raytracer();
}

void XDGUniverse::init_metadata()
{
  xdg_instance_->mesh_manager()->parse_metadata();
}

void XDGUniverse::init_geometry()
{
  int32_t next_cell_id = 0;
  for (const auto& c : model::cells) {
    if (c->id_ > next_cell_id)
      next_cell_id = c->id_;
  }
  cell_idx_offset_ = model::cells.size();
  next_cell_id++;

  for (auto volume : xdg_ptr()->mesh_manager()->volumes()) {
    auto c = std::make_unique<XDGCell>(xdg_ptr(), volume);
    c->id_ =     c->id_ = adjust_geometry_ids_ ? next_cell_id++ : volume;
    c->universe_ = this->id_;

    auto in_map = model::cell_map.find(c->id_);
    if (in_map == model::cell_map.end()) {
      model::cell_map[c->id_] = model::cells.size();
      cell_index_map_[volume] = model::cells.size();
    } else {
      warning(fmt::format("XDG Cell IDs: {}", xdg_ids_for_dim(3)));
      fatal_error(fmt::format(
        "XDG Universe {} contains a cell with ID {}, which "
        "already exists elsewhere in the geometry. Setting auto_geom_ids "
        "to True when defining the XDG Universe may "
        "resolve this issue",
        this->id_, c->id_));
    }

    // --- Materials ---

    // determine volume material assignment
    std::string mat_str;
    if (!xdg_ptr()->mesh_manager()->volume_has_property(volume, xdg::PropertyType::MATERIAL)) {
      fatal_error(fmt::format("XDG Volume {} has no material assignment."), c->id_);
    } else {
      mat_str = xdg_ptr()->mesh_manager()->get_volume_property(volume, xdg::PropertyType::MATERIAL).value;
    }

    to_lower(mat_str);

    assign_material(mat_str, c);

    // ----Temperatures-----

    // if the temperature property is set on the XDG volume, use that
    if (xdg_ptr()->mesh_manager()->volume_has_property(volume, xdg::PropertyType::TEMPERATURE)) {
      // read the temperature from the file (assumed to be in K)
      xdg::Property xdg_temperature = xdg_ptr()->mesh_manager()->get_volume_property(volume, xdg::PropertyType::TEMPERATURE);
      double temperature = std::stod(xdg_temperature.value);
      c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * temperature));
    } else { // otherwise, set the temperature as usual
      const auto& mat = model::materials[model::material_map[c->material_[0]]];
      if (mat->temperature() > 0.0) {
        c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * mat->temperature()));
      } else {
      c->sqrtkT_.push_back(
        std::sqrt(K_BOLTZMANN * settings::temperature_default));
      }
    }
    model::cells.emplace_back(std::move(c));
  }

  // allocate the cell overlap count if necessary
  if (settings::check_overlaps) {
    model::overlap_check_count.resize(model::cells.size(), 0);
  }
  // determine the next surface id
  int32_t next_surf_id = 0;
  for (const auto& s : model::surfaces) {
    if (s->id_ > next_surf_id)
      next_surf_id = s->id_;
  }
  surf_idx_offset_ = model::surfaces.size();
  next_surf_id++;

  for (auto surface : xdg_ptr()->mesh_manager()->surfaces()) {
    auto s = std::make_unique<XDGSurface>(xdg_instance_, surface);
    s->id_ = adjust_geometry_ids_ ? next_surf_id++ : surface;

    s->surf_source_ = contains(settings::source_write_surf_id, s->id_);

    // ----- Boundary Conditions--------
    if (xdg_ptr()->mesh_manager()->surface_has_property(surface, xdg::PropertyType::BOUNDARY_CONDITION)) {
      xdg::Property xdg_bc = xdg_ptr()->mesh_manager()->get_surface_property(surface, xdg::PropertyType::BOUNDARY_CONDITION);
      std::string bc_value = xdg_bc.value;
      to_lower(bc_value);

      if (bc_value.empty() || bc_value == "transmit" ||
          bc_value == "transmission") {
        // set to transmission by default (nullptr)
      } else if (bc_value == "vacuum" || "xdg_boundary") {
        s->bc_ = make_unique<VacuumBC>();
      } else if (bc_value == "reflective" || bc_value == "reflect" ||
                bc_value == "reflecting") {
        s->bc_ = make_unique<ReflectiveBC>();
      } else if (bc_value == "periodic") {
        fatal_error("Periodic boundary condition not supported in XDG geometry.");
      } else {
        fatal_error(fmt::format("Unknown boundary condition \"{}\" specified "
                                "on surface {}",
          bc_value, s->id_));
      }
    }

    // add to global array and map
    auto in_map = model::surface_map.find(s->id_);
    if (in_map == model::surface_map.end()) {
      model::surface_map[s->id_] = model::surfaces.size();
      surface_index_map_[surface] = model::surfaces.size() + 1;
    } else {
      warning(fmt::format("XDG Surface IDs: {}", xdg_ids_for_dim(2)));
      fatal_error(fmt::format("Surface ID {} exists in both Universe {} "
                              "and the CSG geometry.",
        s->id_, this->id_));
    }

    model::surfaces.emplace_back(std::move(s));
  }
}

int32_t XDGUniverse::cell_index(xdg::MeshID volume) const
{
  // return the index of the volume in the XDG instance and then
  // adjust by the offset into the model cells for this XDG universe
  return cell_index_map_.at(volume);
}

int32_t XDGUniverse::surface_index(xdg::MeshID surface) const
{
  // return the index of the surface in the XDG instance and then
  // adjust by the offset into the model cells for this XDG universe
  return surface_index_map_.at(surface);
}

std::string XDGUniverse::xdg_ids_for_dim(int dim) const
{
  // get the ID vector
  std::vector<xdg::MeshID> id_vec;
  if (dim == 3) {
    id_vec = xdg_ptr()->mesh_manager()->volumes();
  } else if (dim == 2) {
    id_vec = xdg_ptr()->mesh_manager()->surfaces();
  } else {
    warning(fmt::format("Invalid dimension {} for XDG ID generation", dim));
    return "";
  }

  // sort the vector of ids
  std::sort(id_vec.begin(), id_vec.end());

  // generate a string representation of the ID range(s)
  std::stringstream out;

  int i = 0;
  int start_id = id_vec[0]; // initialize with first ID
  int stop_id;
  int n_ents = id_vec.size();
  // loop over all cells in the universe
  while (i < n_ents) {

    stop_id = id_vec[i];

    // if the next ID is not in this contiguous set of IDS,
    // figure out how to write the string representing this set
    if (id_vec[i + 1] > stop_id + 1) {

      if (start_id != stop_id) {
        // there are several IDs in a row, print condensed version (i.e. 1-10,
        // 12-20)
        out << start_id << "-" << stop_id;
      } else {
        // only one ID in this contiguous block (i.e. 3, 5, 7, 9)
        out << start_id;
      }
      // insert a comma as long as we aren't in the last ID set
      if (i < n_ents - 1) {
        out << ", ";
      }

      // if we are at the end of a set, set the start ID to the first value
      // in the next set.
      start_id = id_vec[++i];
    }

    i++;
  }

  return out.str();
}

int32_t XDGUniverse::implicit_complement_idx() const
{
  // assume the IPC is at the back of the volumes (for now)
  xdg::MeshID ipc = xdg_ptr()->mesh_manager()->implicit_complement();
  if (ipc == xdg::ID_NONE)
    fatal_error("Implicit complement not found in XDG model");
  return cell_index(ipc);
}

bool XDGUniverse::find_cell(GeometryState& p) const
{
  // if the particle isn't in any of the other XDG
  // cells, place it in the implicit complement
  bool found = Universe::find_cell(p);
  if (!found && model::universe_map[this->id_] != model::root_universe) {
    p.lowest_coord().cell() = implicit_complement_idx();
    found = true;
  }
  return found;
}

void XDGUniverse::to_hdf5(hid_t universes_group) const
{
  // Create a group for this universe.
  auto group = create_group(universes_group, fmt::format("universe {}", id_));

  // Write the geometry representation type.
  write_string(group, "geom_type", "xdg", false);

  // Write other properties of the XDG Universe
  write_string(group, "filename", filename(), false);
  write_attribute(
    group, "auto_geom_ids", static_cast<int>(adjust_geometry_ids_));
  write_attribute(
    group, "auto_mat_ids", static_cast<int>(adjust_material_ids_));

  close_group(group);
}

void XDGUniverse::assign_material(
  std::string& mat_string, std::unique_ptr<XDGCell>& c) const
{
  // material void checks
  if (mat_string == "void" || mat_string == "vacuum" || mat_string == "graveyard") {
    c->material_.push_back(MATERIAL_VOID);
    return;
  }


  bool mat_found_by_name = false;
  // attempt to find a material with a matching name
  to_lower(mat_string);
  for (const auto& m : model::materials) {
    std::string m_name = m->name();
    to_lower(m_name);
    if (mat_string == m_name) {
      // assign the material with that name
      if (!mat_found_by_name) {
        mat_found_by_name = true;
        c->material_.push_back(m->id_);
        // report error if more than one material is found
      } else {
        fatal_error(fmt::format(
          "More than one material found with name '{}'. Please ensure "
          "materials "
          "have unique names if using this property to assign materials.",
          mat_string));
      }
    }
  }

  // if no material was set using a name, assign by id
  if (!mat_found_by_name) {
    bool found_by_id = true;
    try {
      auto id = std::stoi(mat_string);
      if (model::material_map.find(id) == model::material_map.end())
        found_by_id = false;
      c->material_.emplace_back(id);
    } catch (const std::invalid_argument&) {
      found_by_id = false;
    }

    // report failure for failed int conversion or missing material
    if (!found_by_id && !mat_found_by_name)
      fatal_error(
        fmt::format("Material with name/ID '{}' not found for volume (cell) {}",
          mat_string, c->id_));
  }

  if (settings::verbosity >= 10) {
    const auto& m = model::materials[model::material_map.at(c->material_[0])];
    std::stringstream msg;
    msg << "XDG material " << mat_string << " was assigned";
    if (mat_found_by_name) {
      msg << " using material name: " << m->name_;
    } else {
      msg << " using material id: " << m->id_;
    }
    write_message(msg.str(), 10);
  }
}

//==============================================================================
// XDG Cell implementation
//==============================================================================

XDGCell::XDGCell(const std::shared_ptr<xdg::XDG>& xdg_ptr, xdg::MeshID xdg_id)
  : Cell {}, XDGGeometryObject(xdg_ptr, xdg_id)
{
  // TODO: Allow XDG cells to be filled with other geometry
  fill_ = C_NONE;
};

std::pair<double, int32_t> XDGCell::distance(
  Position r, Direction u, int32_t on_surface, GeometryState* p) const
{
  // if we've changed direction or we're not on a surface,
  // reset the history and update last direction
  if (u != p->u_last()) {
    p->u_last() = u;
    p->xdg_prev_elements().clear();
  }
  if (on_surface == 0) {
    p->xdg_prev_elements().clear();
  }

  const auto& univ = model::universes[p->lowest_coord().universe()];

  XDGUniverse* xdg_univ = static_cast<XDGUniverse*>(univ.get());
  if (!xdg_univ)
    fatal_error("XDG call made for particle in a non-XDG universe");

  // initialize to lost particle conditions
  int surf_idx = -1;
  double dist = INFINITY;


  // create the ray
  double pnt[3] = {r.x, r.y, r.z};
  double dir[3] = {u.x, u.y, u.z};
  std::pair<double, xdg::MeshID> result = xdg_ptr()->ray_fire(xdg_id(), pnt, dir, INFTY, xdg::HitOrientation::EXITING, &p->xdg_prev_elements());
  if (result.second > 0) {
    surf_idx = xdg_univ->surface_index(result.second);
    dist = result.first;
    // TODO: reduce complexity of these calls
  } else if (model::cell_map[this->id_] != xdg_univ->implicit_complement_idx() || is_root_universe(xdg_univ->id_)) {
    // surface boundary conditions are ignored for projection plotting, meaning
    // that the particle may move through the graveyard (bounding) volume and
    // into the implicit complement on the other side where no intersection will
    // be found. Treating this as a lost particle is problematic when plotting.
    // Instead, the infinite distance and invalid surface index are returned.
    if (settings::run_mode == RunMode::PLOTTING)
      return {INFTY, -1};

    // the particle should be marked as lost immediately if an intersection
    // isn't found in a volume that is not the implicit complement. In the case
    // that the XDG model is the root universe of the geometry, even a missing
    // intersection in the implicit complement should trigger this condition.
    std::string material_id =
      p->material() == MATERIAL_VOID
        ? "-1 (VOID)"
        : std::to_string(model::materials[p->material()]->id());
    p->mark_as_lost(fmt::format(
      "No intersection found with XDG cell {}, filled with material {}", id_,
      material_id));
  }
  return {dist, surf_idx};
}

bool XDGCell::contains(Position r, Direction u, int32_t on_surface) const
{
  double pnt[3] = {r.x, r.y, r.z};
  xdg::Direction dir {u.x, u.y, u.z};
  return xdg_ptr()->point_in_volume(xdg_id(), pnt, &dir);
}

void XDGCell::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "geom_type", "XDG", false);
}

BoundingBox XDGCell::bounding_box() const
{
  xdg::BoundingBox bbox = xdg_ptr()->mesh_manager()->volume_bounding_box(xdg_id());
  return {{bbox[0], bbox[1], bbox[2]}, {bbox[3], bbox[4], bbox[5]}};
}

//==============================================================================
// XDGSurface implementation
//==============================================================================

XDGSurface::XDGSurface(std::shared_ptr<xdg::XDG> xdg_ptr, int32_t xdg_id)
  : Surface {}, XDGGeometryObject(xdg_ptr, xdg_id)
{}

double XDGSurface::evaluate(Position r) const
{
  return 0.0;
}

double XDGSurface::distance(Position r, Direction u, bool coincident) const
{
  double pnt[3] = {r.x, r.y, r.z};
  double dir[3] = {u.x, u.y, u.z};
  std::pair<double, xdg::MeshID> result = xdg_ptr()->ray_fire(xdg_id(), pnt, dir);
  return result.first < 0.0 ? INFTY : result.first;
}

Direction XDGSurface::normal(Position r) const
{
  xdg::Direction normal = xdg_ptr()->surface_normal(xdg_id(), {r.x, r.y, r.z});
  return {normal[0], normal[1], normal[2]};
}

Direction XDGSurface::reflect(Position r, Direction u, GeometryState* p) const
{
  p->xdg_prev_elements() = {p->xdg_prev_elements().back()};
  double pnt[3] = {r.x, r.y, r.z};
  xdg::Direction normal = xdg_ptr()->surface_normal(xdg_id(), pnt, &p->xdg_prev_elements());
  p->u_last() = u.reflect({normal[0], normal[1], normal[2]});
  return p->u_last();
}

//==============================================================================
// Non-member functions
//==============================================================================

void read_xdg_universes(pugi::xml_node node)
{
  for (pugi::xml_node xdg_node : node.children("xdg_universe")) {
    if (!check_for_node(xdg_node, "type")) {
      fatal_error(fmt::format("XDG type is not specified for universe {}", xdg_node.attribute("id").as_int()));
    }
    std::string type = get_node_value(xdg_node, "type");
    if (type == "surface_mesh") {
      model::universes.push_back(std::make_unique<XDGUniverse>(xdg_node));
    } else if (type == "volume_mesh") {
      model::universes.push_back(std::make_unique<XDGMeshUniverse>(xdg_node));
    } else {
      fatal_error(fmt::format("Invalid XDG universe type: {}", type));
    }
    model::universe_map[model::universes.back()->id_] = model::universes.size() - 1;
  }
}

void check_xdg_root_univ()
{
  return;
  //// TODO: Find an alternative method for this check
  // const auto& ru = model::universes[model::root_universe];
  // if (ru->geom_type() == GeometryType::XDG) {
  //   // if the root universe contains XDG geometry, warn the user
  //   // if it does not contain a graveyard volume
  //   auto dag_univ = dynamic_cast<XDGUniverse*>(ru.get());
  //   if (dag_univ && !dag_univ->has_graveyard()) {
  //     warning(
  //       "No graveyard volume found in the XDG model. "
  //       "This may result in lost particles and rapid simulation failure.");
  //   }
  // }
}

int32_t xdg_next_cell(int32_t surf, int32_t curr_cell, int32_t univ)
{
  auto surfp = dynamic_cast<XDGSurface*>(model::surfaces[surf].get());
  auto cellp = dynamic_cast<XDGCell*>(model::cells[curr_cell].get());
  auto univp = static_cast<XDGUniverse*>(model::universes[univ].get());

  xdg::MeshID next_volume = univp->xdg_ptr()->mesh_manager()->next_volume(cellp->xdg_id(), surfp->xdg_id());
  return univp->cell_index(next_volume);
}

// XDG mesh element universe implementation

XDGMeshUniverse::XDGMeshUniverse(pugi::xml_node node)
{
  if (check_for_node(node, "id")) {
    id_ = std::stoi(get_node_value(node, "id"));
  } else {
    fatal_error("Must specify the ID of the Mesh Universe");
  }

  if (check_for_node(node, "name")) {
    name_ = std::stoi(get_node_value(node, "name"));
  }

  if (check_for_node(node, "mesh")) {
    int32_t mesh_id = std::stoi(get_node_value(node, "mesh"));
    // Ensure the speicifed mesh is present
    if (model::mesh_map.find(mesh_id) == model::mesh_map.end()) {
      fatal_error(fmt::format("Mesh {} could not be found", mesh_id));
    }
    mesh_idx_ = model::mesh_map[mesh_id];
  } else {
    fatal_error(fmt::format("No mesh specified on mesh universe {}", id_));
  }

  create_cells(node);
  mesh()->prepare_for_point_location();
}

int32_t XDGMeshUniverse::match_material(const std::string& material_name) const
{
  int32_t mat_idx {MATERIAL_INVALID};
  if (material_name == "void" || material_name == "vacuum" || material_name == "graveyard") {
    return MATERIAL_VOID;
  }

  // attempt to assign material by name
  mat_idx = get_material_by_name(material_name);
  if (mat_idx != MATERIAL_INVALID) {
    write_message(fmt::format("Assigning material {} by name", material_name), 10);
    return mat_idx;
  }

  // attempt to assign material by ID, treating the string as an integer
  try {
    int32_t mat_id = std::stoi(material_name);
    if (model::material_map.find(mat_id) != model::material_map.end()) {
      return model::material_map.at(mat_id);
    }
  } catch (const std::invalid_argument&) {
    // continue to next step
  }

  fatal_error(fmt::format("Material with name/ID '{}' not found", material_name));
  return MATERIAL_INVALID;
}

void XDGMeshUniverse::create_cells(pugi::xml_node node)
{
  // vector<std::string> cell_fills;
  // if (check_for_node(node, "fills")) {
  //   cell_fills = get_node_array<std::string>(node, "fills");
  // }


  // int n_bins = mesh->n_bins();
  // if (cell_fills.size() != 1 && cell_fills.size() != n_bins) {
  //   fatal_error(
  //     fmt::format("Invalid number of cell fills provided for structured mesh "
  //                 "universe {}. Must be 1 or {}",
  //       id_, n_bins));
  // }

  int32_t next_cell_id = 0;
  for (const auto& c : model::cells) {
    if (c->id_ > next_cell_id)
      next_cell_id = c->id_;
  }
  next_cell_id++;

  cells_.resize(xdg_instance()->mesh_manager()->num_volume_elements(), C_NONE);
  // loop over all volumes in the xdg instance
  for (const auto& volume : xdg_instance()->mesh_manager()->volumes()) {
    const auto& material = xdg_instance()->mesh_manager()->get_volume_property(volume, xdg::PropertyType::MATERIAL);
    std::string mat_str = material.value;
    int32_t mat_idx = match_material(to_lower(mat_str));

    // go over all elements in the volume and assign the material to the element
    for (const auto& element : xdg_instance()->mesh_manager()->get_volume_elements(volume)) {
      element_material_map_[element].push_back(mat_idx);
      // create a cell for the element
      model::cells.push_back(std::make_unique<XDGMeshCell>(mesh_idx_, element));
      auto& c = model::cells.back();
      c->id_ = next_cell_id++;
      model::cell_map[model::cells.back()->id_] = model::cells.size() - 1;

      c->universe_ = id_;
      c->fill_ = C_NONE;
      if (mat_idx != MATERIAL_VOID) {
        c->material_.push_back(model::materials[mat_idx]->id_);
      } else {
        c->material_.push_back(MATERIAL_VOID);
      }
      c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * settings::temperature_default));
      int32_t element_idx = xdg_mesh()->mesh_id_to_bin(element);
      cells_[element_idx] = model::cell_map[c->id_];
    }
  } // end of volume loop

  // create a cell for the exterior of the mesh
  if (check_for_node(node, "background_material")) {
    // get id of outer material
    outer_material_ = match_material(get_node_value(node, "background_material"));
    if (outer_material_ == MATERIAL_INVALID) {
      outer_material_ = MATERIAL_VOID;
    }
  }
  // create a cell for the exterior of the mesh
  model::cells.push_back(std::make_unique<XDGMeshCell>(mesh_idx_, C_NONE));
  auto& c = model::cells.back();
  c->id_ = next_cell_id++;
  model::cell_map[c->id_] = model::cells.size() - 1;
  cells_.push_back(model::cell_map[c->id_]);

  // set exterior cell properties
  c->universe_ = id_;
  c->fill_ = C_NONE;
  if (outer_material_ != MATERIAL_VOID) {
    c->material_.push_back(model::materials[outer_material_]->id_);
  } else {
    c->material_.push_back(MATERIAL_VOID);
  }
  c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * settings::temperature_default));
}

bool XDGMeshUniverse::find_cell(GeometryState& p) const
{
  Position r {p.r_local()};
  bool in_mesh;
  int mesh_bin = mesh()->get_bin(r);

  if (mesh_bin == C_NONE) {
    if (outer_material() == MATERIAL_INVALID)
      return false;
    p.lowest_coord().mesh_cell_index() = C_NONE;
    p.lowest_coord().cell() = cells_[cells_.size() - 1];
    return true;
  }

  p.lowest_coord().mesh_cell_index() = mesh_bin;
  p.lowest_coord().cell() = cells_[mesh_bin];

  return true;
}

void XDGMeshUniverse::next_cell(Particle& p) const
{
  auto& coord = p.lowest_coord();

  int32_t next_mesh_idx = p.boundary().mesh_translation(0);
  int32_t next_cell_idx {C_NONE};
  if (xdg_mesh()->bin_is_valid(next_mesh_idx)) {

    if (p.lowest_coord().mesh_cell_index() == C_NONE) {
      write_message(
        fmt::format(
          "\tParticle {} moving into the mesh. \n\tPosition: {} {} {}", p.id(),
          p.r()[0], p.r()[1], p.r()[2]),
        10);
    }
    write_message(
      fmt::format("\nMoving into mesh cell: {} {} {}",
        p.boundary().mesh_translation(0), p.boundary().mesh_translation(1),
        p.boundary().mesh_translation(2)),
      10);

    next_cell_idx = cells_[next_mesh_idx];
  } else {
    write_message(
      fmt::format("\tParticle {} moving out of the mesh. \n\tPosition: {} {} "
                  "{} \n\tDirection: {} {} {}",
        p.id(), p.r()[0], p.r()[1], p.r()[2], p.u()[0], p.u()[1], p.u()[2]),
      10);
      // treat mesh exit as vacuum boundary for now
      // p.wgt() = 0.0;
      // return;
    next_mesh_idx = C_NONE;
    next_cell_idx = cells_[cells_.size() - 1];
  }

  // reset the lattice_translation for the boundary crossing
  p.boundary().mesh_translation() = {0, 0, 0};
  // update material and temperature
  p.material_last() = p.material();
  p.sqrtkT_last() = p.sqrtkT();
  // set previous bin
  coord.mesh_cell_index() = next_mesh_idx;
  coord.cell() = next_cell_idx;
  coord.mesh_index() = p.boundary().mesh_translation();
  const auto& cell = model::cells.at(next_cell_idx);
  Cell* cell_ptr = model::cells.at(next_cell_idx).get();
  // TODO: Support multiple cell instances
  p.cell_instance() = 0;
  p.material() = cell->material(p.cell_instance());
  p.sqrtkT() = cell->sqrtkT(p.cell_instance());
}

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
  xdg_->prepare_raytracer();
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
  double vol = xdg_->mesh_manager()->element_volume(bin_to_mesh_id(bin));
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

NextMeshCell XDGMesh::distance_to_bin_boundary(GeometryState& g) const
{
  const auto& coord = g.lowest_coord();
  int bin = coord.mesh_cell_index();
  Position r {g.r().x, g.r().y, g.r().z};
  Direction u {g.u().x, g.u().y, g.u().z};
  if (bin == C_NONE || bin == xdg_->mesh_manager()->implicit_complement()) {
    auto ipc = xdg_->mesh_manager()->implicit_complement();
    auto ipc_elem = xdg_->ray_fire(ipc, {r.x, r.y, r.z}, {u.x, u.y, u.z}, INFTY, xdg::HitOrientation::EXITING, &g.xdg_prev_elements());
    if (ipc_elem.second == C_NONE) {
      return {INFTY, -1, {-1, 0, 0}};
    }
    auto new_r = r + u * (ipc_elem.first + TINY_BIT);
    auto next_element = xdg_->find_element({new_r.x, new_r.y, new_r.z});
    return {ipc_elem.first, -1, {mesh_id_to_bin(next_element), 0, 0}};
  }
  auto mesh_id = bin_to_mesh_id(bin);
  auto dist = xdg_->next_element(mesh_id, {r.x, r.y, r.z}, {u.x, u.y, u.z});
  return {dist.second, -1, {mesh_id_to_bin(dist.first), 0, 0}};
}

NextMeshCell XDGMesh::distance_to_bin_boundary(int bin, const Position& r, const Direction& u) const
{
  if (bin == C_NONE || bin == xdg_->mesh_manager()->implicit_complement()) {
    auto ipc = xdg_->mesh_manager()->implicit_complement();
    auto ipc_elem = xdg_->ray_fire(ipc, {r.x, r.y, r.z}, {u.x, u.y, u.z});
    if (ipc_elem.second == C_NONE) {
      return {INFTY, -1, {-1, 0, 0}};
    }
    auto new_r = r + u * (ipc_elem.first + TINY_BIT);
    auto next_element = xdg_->find_element({new_r.x, new_r.y, new_r.z});
    return {ipc_elem.first, -1, {next_element, 0, 0}};
  }
  auto mesh_id = bin_to_mesh_id(bin);
  auto dist = xdg_->next_element(mesh_id, {r.x, r.y, r.z}, {u.x, u.y, u.z});
  return {dist.second, -1, {dist.first, 0, 0}};
}


} // namespace openmc

#else

namespace openmc {

void read_xdg_universes(pugi::xml_node node)
{
  // TODO : updated to xdg universe
  if (check_for_node(node, "xdg_universe")) {
    fatal_error("XDG Universes are present but OpenMC was not configured "
                "with XDG");
  }
};

void check_xdg_root_univ() {};

int32_t next_cell(int32_t surf, int32_t curr_cell, int32_t univ);

} // namespace openmc

#endif // XDG
