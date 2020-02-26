#include "openmc/dagmc.h"

#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/geometry_aux.h"
#include "openmc/material.h"
#include "openmc/string_utils.h"
#include "openmc/settings.h"
#include "openmc/surface.h"

#ifdef DAGMC
#include "uwuw.hpp"
#include "dagmcmetadata.hpp"
#endif
#include <fmt/core.h>

#include <string>
#include <sstream>
#include <algorithm>
#include <fstream>

namespace openmc {

#ifdef DAGMC
const bool dagmc_enabled = true;
#else
const bool dagmc_enabled = false;
#endif

}

#ifdef DAGMC

namespace openmc {

const std::string DAGMC_FILENAME = "dagmc.h5m";

namespace model {

moab::DagMC* DAG;

} // namespace model


void check_dagmc_file() {
  std::string filename = settings::path_input + DAGMC_FILENAME;
  if (!file_exists(filename)) {
    fatal_error("Geometry DAGMC file '" + filename + "' does not exist!");
  }
}

bool get_uwuw_materials_xml(std::string& s) {
  check_dagmc_file();
  UWUW uwuw((settings::path_input + DAGMC_FILENAME).c_str());

  std::stringstream ss;
  bool uwuw_mats_present = false;
  if (uwuw.material_library.size() != 0) {
    uwuw_mats_present = true;
    // write header
    ss << "<?xml version=\"1.0\"?>\n";
    ss << "<materials>\n";
    const auto& mat_lib = uwuw.material_library;
    // write materials
    for (auto mat : mat_lib) { ss << mat.second.openmc("atom"); }
    // write footer
    ss << "</materials>";
    s = ss.str();
  }

  return uwuw_mats_present;
}

bool read_uwuw_materials(pugi::xml_document& doc) {
  std::string s;
  bool found_uwuw_mats = get_uwuw_materials_xml(s);
  if (found_uwuw_mats) {
    pugi::xml_parse_result result = doc.load_string(s.c_str());
  }
  return found_uwuw_mats;
}

bool write_uwuw_materials_xml() {
  std::string s;
  bool found_uwuw_mats = get_uwuw_materials_xml(s);
    // if there is a material library in the file
  if (found_uwuw_mats) {
    // write a material.xml file
    std::ofstream mats_xml("materials.xml");
    mats_xml << s;
    mats_xml.close();
  }

  return found_uwuw_mats;
}

void legacy_assign_material(const std::string& mat_string, DAGCell* c)
{
  bool mat_found_by_name = false;
  // attempt to find a material with a matching name
  for (const auto& m : model::materials) {
    if (mat_string == m->name_) {
      // assign the material with that name
      if (!mat_found_by_name) {
        mat_found_by_name = true;
        c->material_.push_back(m->id_);
      // report error if more than one material is found
      } else {
        fatal_error(fmt::format(
          "More than one material found with name {}. Please ensure materials "
          "have unique names if using this property to assign materials.",
          mat_string));
      }
    }
  }

  // if no material was set using a name, assign by id
  if (!mat_found_by_name) {
    try {
      auto id = std::stoi(mat_string);
      c->material_.emplace_back(id);
    } catch (const std::invalid_argument&) {
      fatal_error(fmt::format(
        "No material {} found for volume (cell) {}", mat_string, c->id_));
    }
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

void load_dagmc_geometry()
{
  check_dagmc_file();

  if (!model::DAG) {
    model::DAG = new moab::DagMC();
  }

  std::string filename = settings::path_input + DAGMC_FILENAME;
  // --- Materials ---

  // create uwuw instance
  UWUW uwuw(filename.c_str());

  // check for uwuw material definitions
  bool using_uwuw = !uwuw.material_library.empty();

  // notify user if UWUW materials are going to be used
  if (using_uwuw) {
    write_message("Found UWUW Materials in the DAGMC geometry file.", 6);
  }

  int32_t dagmc_univ_id = 0; // universe is always 0 for DAGMC runs

  // load the DAGMC geometry
  moab::ErrorCode rval = model::DAG->load_file(filename.c_str());
  MB_CHK_ERR_CONT(rval);

  // initialize acceleration data structures
  rval = model::DAG->init_OBBTree();
  MB_CHK_ERR_CONT(rval);

  // parse model metadata
  dagmcMetaData DMD(model::DAG);
  if (using_uwuw) {
    DMD.load_property_data();
  }
  std::vector<std::string> keywords {"temp", "mat", "density", "boundary"};
  std::map<std::string, std::string> dum;
  std::string delimiters = ":/";
  rval = model::DAG->parse_properties(keywords, dum, delimiters.c_str());
  MB_CHK_ERR_CONT(rval);

  // --- Cells (Volumes) ---

  // initialize cell objects
  int n_cells = model::DAG->num_entities(3);
  moab::EntityHandle graveyard = 0;
  for (int i = 0; i < n_cells; i++) {
    moab::EntityHandle vol_handle = model::DAG->entity_by_index(3, i+1);

    // set cell ids using global IDs
    DAGCell* c = new DAGCell();
    c->dag_index_ = i+1;
    c->id_ = model::DAG->id_by_index(3, c->dag_index_);
    c->dagmc_ptr_ = model::DAG;
    c->universe_ = dagmc_univ_id; // set to zero for now
    c->fill_ = C_NONE; // no fill, single universe

    model::cells.emplace_back(c);
    model::cell_map[c->id_] = i;

    // Populate the Universe vector and dict
    auto it = model::universe_map.find(dagmc_univ_id);
    if (it == model::universe_map.end()) {
      model::universes.push_back(std::make_unique<Universe>());
      model::universes.back()->id_ = dagmc_univ_id;
      model::universes.back()->cells_.push_back(i);
      model::universe_map[dagmc_univ_id] = model::universes.size() - 1;
    } else {
      model::universes[it->second]->cells_.push_back(i);
    }

    // MATERIALS

    if (model::DAG->is_implicit_complement(vol_handle)) {
      if (model::DAG->has_prop(vol_handle, "mat")) {
        // if the implicit complement has been assigned a material, use it
        if (using_uwuw) {
          std::string comp_mat = DMD.volume_material_property_data_eh[vol_handle];
          // Note: material numbers are set by UWUW
          int mat_number = uwuw.material_library[comp_mat].metadata["mat_number"].asInt();
          c->material_.push_back(mat_number);
        } else {
          std::string mat_value;
          rval= model::DAG->prop_value(vol_handle, "mat", mat_value);
          MB_CHK_ERR_CONT(rval);
          // remove _comp
          std::string _comp = "_comp";
          size_t _comp_pos = mat_value.find(_comp);
          if (_comp_pos != std::string::npos) { mat_value.erase(_comp_pos, _comp.length()); }
          // assign IC material by id
          legacy_assign_material(mat_value, c);
        }
      } else {
        // if no material is found, the implicit complement is void
        c->material_.push_back(MATERIAL_VOID);
      }
      continue;
    }

    // determine volume material assignment
    std::string mat_value;
    if (model::DAG->has_prop(vol_handle, "mat")) {
      rval = model::DAG->prop_value(vol_handle, "mat", mat_value);
      MB_CHK_ERR_CONT(rval);
    } else {
      fatal_error(fmt::format("Volume {} has no material assignment.", c->id_));
    }

    std::string cmp_str = mat_value;
    to_lower(cmp_str);

    if (cmp_str == "graveyard") {
      graveyard = vol_handle;
    }

    // material void checks
    if (cmp_str == "void" || cmp_str == "vacuum" || cmp_str == "graveyard") {
      c->material_.push_back(MATERIAL_VOID);
    } else {
      if (using_uwuw) {
        // lookup material in uwuw if the were present
        std::string uwuw_mat = DMD.volume_material_property_data_eh[vol_handle];
        if (uwuw.material_library.count(uwuw_mat) != 0) {
          // Note: material numbers are set by UWUW
          int mat_number = uwuw.material_library[uwuw_mat].metadata["mat_number"].asInt();
          c->material_.push_back(mat_number);
        } else {
          fatal_error(fmt::format("Material with value {} not found in the "
            "UWUW material library", mat_value));
        }
      } else {
        legacy_assign_material(mat_value, c);
      }
    }

    // check for temperature assignment
    std::string temp_value;

    // no temperature if void
    if (c->material_[0] == MATERIAL_VOID) continue;

    // assign cell temperature
    const auto& mat = model::materials[model::material_map.at(c->material_[0])];
    if (model::DAG->has_prop(vol_handle, "temp")) {
      rval = model::DAG->prop_value(vol_handle, "temp", temp_value);
      MB_CHK_ERR_CONT(rval);
      double temp = std::stod(temp_value);
      c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * temp));
    } else if (mat->temperature_ > 0.0) {
      c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * mat->temperature_));
    } else {
      c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * settings::temperature_default));
    }

  }

  // allocate the cell overlap count if necessary
  if (settings::check_overlaps) {
    model::overlap_check_count.resize(model::cells.size(), 0);
  }

  if (!graveyard) {
    warning("No graveyard volume found in the DagMC model."
            "This may result in lost particles and rapid simulation failure.");
  }

  // --- Surfaces ---

  // initialize surface objects
  int n_surfaces = model::DAG->num_entities(2);

  for (int i = 0; i < n_surfaces; i++) {
    moab::EntityHandle surf_handle = model::DAG->entity_by_index(2, i+1);

    // set cell ids using global IDs
    DAGSurface* s = new DAGSurface();
    s->dag_index_ = i+1;
    s->id_ = model::DAG->id_by_index(2, s->dag_index_);
    s->dagmc_ptr_ = model::DAG;

    // set BCs
    std::string bc_value;
    if (model::DAG->has_prop(surf_handle, "boundary")) {
      rval = model::DAG->prop_value(surf_handle, "boundary", bc_value);
      MB_CHK_ERR_CONT(rval);
      to_lower(bc_value);

      if (bc_value == "transmit" || bc_value == "transmission") {
        s->bc_ = Surface::BoundaryType::TRANSMIT;
      } else if (bc_value == "vacuum") {
        s->bc_ = Surface::BoundaryType::VACUUM;
      } else if (bc_value == "reflective" || bc_value == "reflect" || bc_value == "reflecting") {
        s->bc_ = Surface::BoundaryType::REFLECT;
      } else if (bc_value == "periodic") {
        fatal_error("Periodic boundary condition not supported in DAGMC.");
      } else {
        fatal_error(fmt::format("Unknown boundary condition \"{}\" specified "
          "on surface {}", bc_value, s->id_));
      }
    } else {
      // if no condition is found, set to transmit
      s->bc_ = Surface::BoundaryType::TRANSMIT;
    }

    // graveyard check
    moab::Range parent_vols;
    rval = model::DAG->moab_instance()->get_parent_meshsets(surf_handle, parent_vols);
    MB_CHK_ERR_CONT(rval);

    // if this surface belongs to the graveyard
    if (graveyard && parent_vols.find(graveyard) != parent_vols.end()) {
      // set BC to vacuum
      s->bc_ = Surface::BoundaryType::VACUUM;
    }

    // add to global array and map
    model::surfaces.emplace_back(s);
    model::surface_map[s->id_] = i;
  }

  return;
}

void read_geometry_dagmc()
{
  // Check if dagmc.h5m exists
  check_dagmc_file();
  write_message("Reading DAGMC geometry...", 5);
  load_dagmc_geometry();

  model::root_universe = find_root_universe();
}

void free_memory_dagmc()
{
  delete model::DAG;
}


}
#endif
