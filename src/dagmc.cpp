#include "openmc/dagmc.h"

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

#ifdef UWUW_HPP
#include "uwuw.hpp"
#endif
#include <fmt/core.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

namespace openmc {

#ifdef DAGMC
const bool DAGMC_ENABLED = true;
#else
const bool DAGMC_ENABLED = false;
#endif

#ifdef UWUW_HPP
const bool UWUW_ENABLED = true;
#else
const bool UWUW_ENABLED = false;
#endif

} // namespace openmc

#ifdef DAGMC

namespace openmc {

//==============================================================================
// DAGMC Universe implementation
//==============================================================================

DAGUniverse::DAGUniverse(pugi::xml_node node)
{
  if (check_for_node(node, "id")) {
    id_ = std::stoi(get_node_value(node, "id"));
  } else {
    fatal_error("Must specify the id of the DAGMC universe");
  }

  if (check_for_node(node, "filename")) {
    filename_ = get_node_value(node, "filename");
    if (!starts_with(filename_, "/")) {
      filename_ = dir_name(settings::path_input) + filename_;
    }
  } else {
    fatal_error("Must specify a file for the DAGMC universe");
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

DAGUniverse::DAGUniverse(
  const std::string& filename, bool auto_geom_ids, bool auto_mat_ids)
  : filename_(filename), adjust_geometry_ids_(auto_geom_ids),
    adjust_material_ids_(auto_mat_ids)
{
  set_id();
  initialize();
}

DAGUniverse::DAGUniverse(std::shared_ptr<moab::DagMC> dagmc_ptr,
  const std::string& filename, bool auto_geom_ids, bool auto_mat_ids)
  : dagmc_instance_(dagmc_ptr), filename_(filename),
    adjust_geometry_ids_(auto_geom_ids), adjust_material_ids_(auto_mat_ids)
{
  set_id();
  init_metadata();
  init_geometry();
}

void DAGUniverse::set_id()
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

void DAGUniverse::initialize()
{
  geom_type() = GeometryType::DAG;

  init_dagmc();

  init_metadata();

  init_geometry();
}

void DAGUniverse::init_dagmc()
{

  // create a new DAGMC instance
  dagmc_instance_ = std::make_shared<moab::DagMC>();

  // load the DAGMC geometry
  if (!file_exists(filename_)) {
    fatal_error("Geometry DAGMC file '" + filename_ + "' does not exist!");
  }
  moab::ErrorCode rval = dagmc_instance_->load_file(filename_.c_str());
  MB_CHK_ERR_CONT(rval);

  // initialize acceleration data structures
  rval = dagmc_instance_->init_OBBTree();
  MB_CHK_ERR_CONT(rval);
}

void DAGUniverse::init_metadata()
{
  // parse model metadata
  dmd_ptr =
    std::make_unique<dagmcMetaData>(dagmc_instance_.get(), false, false);
  dmd_ptr->load_property_data();

  std::vector<std::string> keywords {"temp"};
  std::map<std::string, std::string> dum;
  std::string delimiters = ":/";
  moab::ErrorCode rval;
  rval = dagmc_instance_->parse_properties(keywords, dum, delimiters.c_str());
  MB_CHK_ERR_CONT(rval);
}

void DAGUniverse::init_geometry()
{
  moab::ErrorCode rval;

  // determine the next cell id
  int32_t next_cell_id = 0;
  for (const auto& c : model::cells) {
    if (c->id_ > next_cell_id)
      next_cell_id = c->id_;
  }
  cell_idx_offset_ = model::cells.size();
  next_cell_id++;

  // initialize cell objects
  int n_cells = dagmc_instance_->num_entities(3);
  moab::EntityHandle graveyard = 0;
  for (int i = 0; i < n_cells; i++) {
    moab::EntityHandle vol_handle = dagmc_instance_->entity_by_index(3, i + 1);

    // set cell ids using global IDs
    auto c = std::make_unique<DAGCell>(dagmc_instance_, i + 1);
    c->id_ = adjust_geometry_ids_
               ? next_cell_id++
               : dagmc_instance_->id_by_index(3, c->dag_index());
    c->universe_ = this->id_;
    c->fill_ = C_NONE; // no fill, single universe

    auto in_map = model::cell_map.find(c->id_);
    if (in_map == model::cell_map.end()) {
      model::cell_map[c->id_] = model::cells.size();
    } else {
      warning(fmt::format("DAGMC Cell IDs: {}", dagmc_ids_for_dim(3)));
      fatal_error(fmt::format(
        "DAGMC Universe {} contains a cell with ID {}, which "
        "already exists elsewhere in the geometry. Setting auto_geom_ids "
        "to True when initiating the DAGMC Universe may "
        "resolve this issue",
        this->id_, c->id_));
    }

    // --- Materials ---

    // determine volume material assignment
    std::string mat_str = dmd_ptr->get_volume_property("material", vol_handle);

    if (mat_str.empty()) {
      fatal_error(fmt::format("Volume {} has no material assignment.", c->id_));
    }

    to_lower(mat_str);

    if (mat_str == "graveyard") {
      graveyard = vol_handle;
    }

    // material void checks
    if (mat_str == "void" || mat_str == "vacuum" || mat_str == "graveyard") {
      c->material_.push_back(MATERIAL_VOID);
    } else {
      if (uses_uwuw()) {
        uwuw_assign_material(vol_handle, c);
      } else {
        legacy_assign_material(mat_str, c);
      }
    }

    // check for temperature assignment
    std::string temp_value;

    // no temperature if void
    if (c->material_[0] == MATERIAL_VOID) {
      model::cells.emplace_back(std::move(c));
      continue;
    }

    // assign cell temperature
    const auto& mat = model::materials[model::material_map.at(c->material_[0])];
    if (dagmc_instance_->has_prop(vol_handle, "temp")) {
      rval = dagmc_instance_->prop_value(vol_handle, "temp", temp_value);
      MB_CHK_ERR_CONT(rval);
      double temp = std::stod(temp_value);
      c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * temp));
    } else if (mat->temperature() > 0.0) {
      c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * mat->temperature()));
    } else {
      c->sqrtkT_.push_back(
        std::sqrt(K_BOLTZMANN * settings::temperature_default));
    }

    model::cells.emplace_back(std::move(c));
  }

  // allocate the cell overlap count if necessary
  if (settings::check_overlaps) {
    model::overlap_check_count.resize(model::cells.size(), 0);
  }

  has_graveyard_ = graveyard;

  // determine the next surface id
  int32_t next_surf_id = 0;
  for (const auto& s : model::surfaces) {
    if (s->id_ > next_surf_id)
      next_surf_id = s->id_;
  }
  surf_idx_offset_ = model::surfaces.size();
  next_surf_id++;

  // initialize surface objects
  int n_surfaces = dagmc_instance_->num_entities(2);
  for (int i = 0; i < n_surfaces; i++) {
    moab::EntityHandle surf_handle = dagmc_instance_->entity_by_index(2, i + 1);

    // set cell ids using global IDs
    auto s = std::make_unique<DAGSurface>(dagmc_instance_, i + 1);
    s->id_ = adjust_geometry_ids_ ? next_surf_id++
                                  : dagmc_instance_->id_by_index(2, i + 1);

    // set surface source attribute if needed
    if (contains(settings::source_write_surf_id, s->id_) ||
        settings::source_write_surf_id.empty()) {
      s->surf_source_ = true;
    }

    // set BCs
    std::string bc_value =
      dmd_ptr->get_surface_property("boundary", surf_handle);
    to_lower(bc_value);
    if (bc_value.empty() || bc_value == "transmit" ||
        bc_value == "transmission") {
      // set to transmission by default (nullptr)
    } else if (bc_value == "vacuum") {
      s->bc_ = make_unique<VacuumBC>();
    } else if (bc_value == "reflective" || bc_value == "reflect" ||
               bc_value == "reflecting") {
      s->bc_ = make_unique<ReflectiveBC>();
    } else if (bc_value == "periodic") {
      fatal_error("Periodic boundary condition not supported in DAGMC.");
    } else {
      fatal_error(fmt::format("Unknown boundary condition \"{}\" specified "
                              "on surface {}",
        bc_value, s->id_));
    }

    // graveyard check
    moab::Range parent_vols;
    rval = dagmc_instance_->moab_instance()->get_parent_meshsets(
      surf_handle, parent_vols);
    MB_CHK_ERR_CONT(rval);

    // if this surface belongs to the graveyard
    if (graveyard && parent_vols.find(graveyard) != parent_vols.end()) {
      // set graveyard surface BC's to vacuum
      s->bc_ = make_unique<VacuumBC>();
    }

    // add to global array and map

    auto in_map = model::surface_map.find(s->id_);
    if (in_map == model::surface_map.end()) {
      model::surface_map[s->id_] = model::surfaces.size();
    } else {
      warning(fmt::format("DAGMC Surface IDs: {}", dagmc_ids_for_dim(2)));
      fatal_error(fmt::format("Surface ID {} exists in both Universe {} "
                              "and the CSG geometry.",
        s->id_, this->id_));
    }

    model::surfaces.emplace_back(std::move(s));
  } // end surface loop
}

int32_t DAGUniverse::cell_index(moab::EntityHandle vol) const
{
  // return the index of the volume in the DAGMC instance and then
  // adjust by the offset into the model cells for this DAGMC universe
  return dagmc_ptr()->index_by_handle(vol) + cell_idx_offset_;
}

int32_t DAGUniverse::surface_index(moab::EntityHandle surf) const
{
  // return the index of the surface in the DAGMC instance and then
  // adjust by the offset into the model cells for this DAGMC universe
  return dagmc_ptr()->index_by_handle(surf) + surf_idx_offset_;
}

std::string DAGUniverse::dagmc_ids_for_dim(int dim) const
{
  // generate a vector of ids
  std::vector<int> id_vec;
  int n_ents = dagmc_instance_->num_entities(dim);
  for (int i = 1; i <= n_ents; i++) {
    id_vec.push_back(dagmc_instance_->id_by_index(dim, i));
  }

  // sort the vector of ids
  std::sort(id_vec.begin(), id_vec.end());

  // generate a string representation of the ID range(s)
  std::stringstream out;

  int i = 0;
  int start_id = id_vec[0]; // initialize with first ID
  int stop_id;
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

int32_t DAGUniverse::implicit_complement_idx() const
{
  moab::EntityHandle ic;
  moab::ErrorCode rval =
    dagmc_instance_->geom_tool()->get_implicit_complement(ic);
  MB_CHK_SET_ERR_CONT(rval, "Failed to get implicit complement");
  // off-by-one: DAGMC indices start at one
  return cell_idx_offset_ + dagmc_instance_->index_by_handle(ic) - 1;
}

bool DAGUniverse::find_cell(GeometryState& p) const
{
  // if the particle isn't in any of the other DagMC
  // cells, place it in the implicit complement
  bool found = Universe::find_cell(p);
  if (!found && model::universe_map[this->id_] != model::root_universe) {
    p.lowest_coord().cell = implicit_complement_idx();
    found = true;
  }
  return found;
}

void DAGUniverse::to_hdf5(hid_t universes_group) const
{
  // Create a group for this universe.
  auto group = create_group(universes_group, fmt::format("universe {}", id_));

  // Write the geometry representation type.
  write_string(group, "geom_type", "dagmc", false);

  // Write other properties of the DAGMC Universe
  write_string(group, "filename", filename_, false);
  write_attribute(
    group, "auto_geom_ids", static_cast<int>(adjust_geometry_ids_));
  write_attribute(
    group, "auto_mat_ids", static_cast<int>(adjust_material_ids_));

  close_group(group);
}

bool DAGUniverse::uses_uwuw() const
{
#ifdef UWUW_HPP
  return uwuw_ && !uwuw_->material_library.empty();
#else
  return false;
#endif // UWUW_HPP
}

std::string DAGUniverse::get_uwuw_materials_xml() const
{
#ifdef UWUW_HPP
  if (!uses_uwuw()) {
    throw std::runtime_error("This DAGMC Universe does not use UWUW materials");
  }

  std::stringstream ss;
  // write header
  ss << "<?xml version=\"1.0\"?>\n";
  ss << "<materials>\n";
  const auto& mat_lib = uwuw_->material_library;
  // write materials
  for (auto mat : mat_lib) {
    ss << mat.second->openmc("atom");
  }
  // write footer
  ss << "</materials>";

  return ss.str();
#else
  fatal_error("DAGMC was not configured with UWUW.");
#endif // UWUW_HPP
}

void DAGUniverse::write_uwuw_materials_xml(const std::string& outfile) const
{
#ifdef UWUW_HPP
  if (!uses_uwuw()) {
    throw std::runtime_error(
      "This DAGMC universe does not use UWUW materials.");
  }

  std::string xml_str = get_uwuw_materials_xml();
  // if there is a material library in the file
  std::ofstream mats_xml(outfile);
  mats_xml << xml_str;
  mats_xml.close();
#else
  fatal_error("DAGMC was not configured with UWUW.");
#endif // UWUW_HPP
}

void DAGUniverse::legacy_assign_material(
  std::string mat_string, std::unique_ptr<DAGCell>& c) const
{
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
    if (!found_by_id)
      fatal_error(
        fmt::format("Material with name/ID '{}' not found for volume (cell) {}",
          mat_string, c->id_));
  }

  if (settings::verbosity >= 10) {
    const auto& m = model::materials[model::material_map.at(c->material_[0])];
    std::stringstream msg;
    msg << "DAGMC material " << mat_string << " was assigned";
    if (mat_found_by_name) {
      msg << " using material name: " << m->name_;
    } else {
      msg << " using material id: " << m->id_;
    }
    write_message(msg.str(), 10);
  }
}

void DAGUniverse::read_uwuw_materials()
{
#ifdef UWUW_HPP
  // If no filename was provided, don't read UWUW materials
  if (filename_ == "")
    return;

  uwuw_ = std::make_shared<UWUW>(filename_.c_str());

  if (!uses_uwuw())
    return;

  // Notify user if UWUW materials are going to be used
  write_message("Found UWUW Materials in the DAGMC geometry file.", 6);

  // if we're using automatic IDs, update the UWUW material metadata
  if (adjust_material_ids_) {
    int32_t next_material_id = 0;
    for (const auto& m : model::materials) {
      next_material_id = std::max(m->id_, next_material_id);
    }
    next_material_id++;

    for (auto& mat : uwuw_->material_library) {
      mat.second->metadata["mat_number"] = next_material_id++;
    }
  }

  std::string mat_xml_string = get_uwuw_materials_xml();

  // create a pugi XML document from this string
  pugi::xml_document doc;
  auto result = doc.load_string(mat_xml_string.c_str());
  if (!result) {
    fatal_error("Error processing XML created using DAGMC UWUW materials.");
  }
  pugi::xml_node root = doc.document_element();
  for (pugi::xml_node material_node : root.children("material")) {
    model::materials.push_back(std::make_unique<Material>(material_node));
  }
#else
  fatal_error("DAGMC was not configured with UWUW.");
#endif // UWUW_HPP
}

void DAGUniverse::uwuw_assign_material(
  moab::EntityHandle vol_handle, std::unique_ptr<DAGCell>& c)
{
#ifdef UWUW_HPP
  // read materials from uwuw material file
  read_uwuw_materials();

  // lookup material in uwuw if present
  std::string uwuw_mat = dmd_ptr->volume_material_property_data_eh[vol_handle];
  if (uwuw_->material_library.count(uwuw_mat) != 0) {
    // Note: material numbers are set by UWUW
    int mat_number = uwuw_->material_library.get_material(uwuw_mat)
                       .metadata["mat_number"]
                       .asInt();
    c->material_.push_back(mat_number);
  } else {
    fatal_error(fmt::format("Material with value '{}' not found in the "
                            "UWUW material library",
                            uwuw_mat));  // Replaced mat_str with uwuw_mat
  }
#else
  fatal_error("DAGMC was not configured with UWUW.");
#endif // UWUW_HPP
}
//==============================================================================
// DAGMC Cell implementation
//==============================================================================

DAGCell::DAGCell(std::shared_ptr<moab::DagMC> dag_ptr, int32_t dag_idx)
  : Cell {}, dagmc_ptr_(dag_ptr), dag_index_(dag_idx)
{
  geom_type_ = GeometryType::DAG;
};

std::pair<double, int32_t> DAGCell::distance(
  Position r, Direction u, int32_t on_surface, GeometryState* p) const
{
  // if we've changed direction or we're not on a surface,
  // reset the history and update last direction
  if (u != p->last_dir()) {
    p->last_dir() = u;
    p->history().reset();
  }
  if (on_surface == 0) {
    p->history().reset();
  }

  const auto& univ = model::universes[p->lowest_coord().universe];

  DAGUniverse* dag_univ = static_cast<DAGUniverse*>(univ.get());
  if (!dag_univ)
    fatal_error("DAGMC call made for particle in a non-DAGMC universe");

  // initialize to lost particle conditions
  int surf_idx = -1;
  double dist = INFINITY;

  moab::EntityHandle vol = dagmc_ptr_->entity_by_index(3, dag_index_);
  moab::EntityHandle hit_surf;

  // create the ray
  double pnt[3] = {r.x, r.y, r.z};
  double dir[3] = {u.x, u.y, u.z};
  MB_CHK_ERR_CONT(
    dagmc_ptr_->ray_fire(vol, pnt, dir, hit_surf, dist, &p->history()));
  if (hit_surf != 0) {
    surf_idx =
      dag_univ->surf_idx_offset_ + dagmc_ptr_->index_by_handle(hit_surf);
  } else if (!dagmc_ptr_->is_implicit_complement(vol) ||
             is_root_universe(dag_univ->id_)) {
    // surface boundary conditions are ignored for projection plotting, meaning
    // that the particle may move through the graveyard (bounding) volume and
    // into the implicit complement on the other side where no intersection will
    // be found. Treating this as a lost particle is problematic when plotting.
    // Instead, the infinite distance and invalid surface index are returned.
    if (settings::run_mode == RunMode::PLOTTING)
      return {INFTY, -1};

    // the particle should be marked as lost immediately if an intersection
    // isn't found in a volume that is not the implicit complement. In the case
    // that the DAGMC model is the root universe of the geometry, even a missing
    // intersection in the implicit complement should trigger this condition.
    std::string material_id =
      p->material() == MATERIAL_VOID
        ? "-1 (VOID)"
        : std::to_string(model::materials[p->material()]->id());
    p->mark_as_lost(fmt::format(
      "No intersection found with DAGMC cell {}, filled with material {}", id_,
      material_id));
  }

  return {dist, surf_idx};
}

bool DAGCell::contains(Position r, Direction u, int32_t on_surface) const
{
  moab::ErrorCode rval;
  moab::EntityHandle vol = dagmc_ptr_->entity_by_index(3, dag_index_);

  int result = 0;
  double pnt[3] = {r.x, r.y, r.z};
  double dir[3] = {u.x, u.y, u.z};
  rval = dagmc_ptr_->point_in_volume(vol, pnt, result, dir);
  MB_CHK_ERR_CONT(rval);
  return result;
}

moab::EntityHandle DAGCell::mesh_handle() const
{
  return dagmc_ptr()->entity_by_index(3, dag_index());
}

void DAGCell::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "geom_type", "dagmc", false);
}

BoundingBox DAGCell::bounding_box() const
{
  moab::ErrorCode rval;
  moab::EntityHandle vol = dagmc_ptr_->entity_by_index(3, dag_index_);
  double min[3], max[3];
  rval = dagmc_ptr_->getobb(vol, min, max);
  MB_CHK_ERR_CONT(rval);
  return {min[0], max[0], min[1], max[1], min[2], max[2]};
}

//==============================================================================
// DAGSurface implementation
//==============================================================================

DAGSurface::DAGSurface(std::shared_ptr<moab::DagMC> dag_ptr, int32_t dag_idx)
  : Surface {}, dagmc_ptr_(dag_ptr), dag_index_(dag_idx)
{
  geom_type_ = GeometryType::DAG;
} // empty constructor

moab::EntityHandle DAGSurface::mesh_handle() const
{
  return dagmc_ptr()->entity_by_index(2, dag_index());
}

double DAGSurface::evaluate(Position r) const
{
  return 0.0;
}

double DAGSurface::distance(Position r, Direction u, bool coincident) const
{
  moab::ErrorCode rval;
  moab::EntityHandle surf = dagmc_ptr_->entity_by_index(2, dag_index_);
  moab::EntityHandle hit_surf;
  double dist;
  double pnt[3] = {r.x, r.y, r.z};
  double dir[3] = {u.x, u.y, u.z};
  rval = dagmc_ptr_->ray_fire(surf, pnt, dir, hit_surf, dist, NULL, 0, 0);
  MB_CHK_ERR_CONT(rval);
  if (dist < 0.0)
    dist = INFTY;
  return dist;
}

Direction DAGSurface::normal(Position r) const
{
  moab::ErrorCode rval;
  moab::EntityHandle surf = dagmc_ptr_->entity_by_index(2, dag_index_);
  double pnt[3] = {r.x, r.y, r.z};
  double dir[3];
  rval = dagmc_ptr_->get_angle(surf, pnt, dir);
  MB_CHK_ERR_CONT(rval);
  return dir;
}

Direction DAGSurface::reflect(Position r, Direction u, GeometryState* p) const
{
  Expects(p);
  p->history().reset_to_last_intersection();
  moab::ErrorCode rval;
  moab::EntityHandle surf = dagmc_ptr_->entity_by_index(2, dag_index_);
  double pnt[3] = {r.x, r.y, r.z};
  double dir[3];
  rval = dagmc_ptr_->get_angle(surf, pnt, dir, &p->history());
  MB_CHK_ERR_CONT(rval);
  p->last_dir() = u.reflect(dir);
  return p->last_dir();
}

//==============================================================================
// Non-member functions
//==============================================================================

void read_dagmc_universes(pugi::xml_node node)
{
  for (pugi::xml_node dag_node : node.children("dagmc_universe")) {
    model::universes.push_back(std::make_unique<DAGUniverse>(dag_node));
    model::universe_map[model::universes.back()->id_] =
      model::universes.size() - 1;
  }
}

void check_dagmc_root_univ()
{
  const auto& ru = model::universes[model::root_universe];
  if (ru->geom_type() == GeometryType::DAG) {
    // if the root universe contains DAGMC geometry, warn the user
    // if it does not contain a graveyard volume
    auto dag_univ = dynamic_cast<DAGUniverse*>(ru.get());
    if (dag_univ && !dag_univ->has_graveyard()) {
      warning(
        "No graveyard volume found in the DagMC model. "
        "This may result in lost particles and rapid simulation failure.");
    }
  }
}

int32_t next_cell(int32_t surf, int32_t curr_cell, int32_t univ)
{
  auto surfp = dynamic_cast<DAGSurface*>(model::surfaces[surf - 1].get());
  auto cellp = dynamic_cast<DAGCell*>(model::cells[curr_cell].get());
  auto univp = static_cast<DAGUniverse*>(model::universes[univ].get());

  moab::EntityHandle surf_handle = surfp->mesh_handle();
  moab::EntityHandle curr_vol = cellp->mesh_handle();

  moab::EntityHandle new_vol;
  moab::ErrorCode rval =
    cellp->dagmc_ptr()->next_vol(surf_handle, curr_vol, new_vol);
  if (rval != moab::MB_SUCCESS)
    return -1;

  return univp->cell_index(new_vol);
}

} // namespace openmc

#else

namespace openmc {

void read_dagmc_universes(pugi::xml_node node)
{
  if (check_for_node(node, "dagmc_universe")) {
    fatal_error("DAGMC Universes are present but OpenMC was not configured "
                "with DAGMC");
  }
};

void check_dagmc_root_univ() {};

int32_t next_cell(int32_t surf, int32_t curr_cell, int32_t univ);

} // namespace openmc

#endif // DAGMC
