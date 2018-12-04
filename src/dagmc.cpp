#include "openmc/dagmc.h"

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/string_utils.h"
#include "openmc/settings.h"
#include "openmc/geometry.h"

#include <string>
#include <sstream>
#include <algorithm>
#include <fstream>

#ifdef DAGMC

#include "uwuw.hpp"
#include "dagmcmetadata.hpp"

#define DAGMC_FILENAME "dagmc.h5m"

namespace openmc {

namespace model {

moab::DagMC* DAG;

} // namespace model

bool write_materials_xml(UWUW uwuw) {
  // if there is a material library in the file,
  // write a material.xml file
  std::ofstream mats_xml("materials.xml");
  // write header
  mats_xml << "<?xml version=\"1.0\"?>\n";
  mats_xml << "<materials>\n";
  std::map<std::string, pyne::Material> ml = uwuw.material_library;
  std::map<std::string, pyne::Material>::iterator it;
  // write materials
  for (it = ml.begin(); it != ml.end(); it++) { mats_xml << it->second.openmc("atom");  }
  // write footer
  mats_xml << "</materials>";
  mats_xml.close();
}

void load_dagmc_geometry()
{
  if (!model::DAG) {
    model::DAG = new moab::DagMC();
  }

  // create uwuw instance
  UWUW uwuw(DAGMC_FILENAME);

  // check for uwuw material definitions
  bool using_uwuw = (uwuw.material_library.size() == 0) ? false : true;

  if (using_uwuw) {
    std::cout << "Found UWUW Materials in the DAGMC geometry file." << std::endl;
    // don't overwrite an existing materials.xml file
    if (file_exists("materials.xml")) {
        std::stringstream err_msg;
        err_msg << "A materials.xml file is present along with UWUW material definitions";
        fatal_error(err_msg.str());
      }
    write_materials_xml(uwuw);
  }

  int32_t dagmc_univ_id = 0; // universe is always 0 for DAGMC runs

  moab::ErrorCode rval = model::DAG->load_file(DAGMC_FILENAME);
  MB_CHK_ERR_CONT(rval);

  rval = model::DAG->init_OBBTree();
  MB_CHK_ERR_CONT(rval);

  dagmcMetaData DMD(model::DAG);
  if (using_uwuw) {
    DMD.load_property_data();
  }
  
  std::vector<std::string> keywords;
  keywords.push_back("mat");
  keywords.push_back("density");
  keywords.push_back("boundary");
  std::map<std::string, std::string> dum;
  rval = model::DAG->parse_properties(keywords, dum, ":/");
  
  // initialize cell objects
  model::n_cells = model::DAG->num_entities(3);

  for (int i = 0; i < model::n_cells; i++) {
    moab::EntityHandle vol_handle = model::DAG->entity_by_index(3, i+1);

    // set cell ids using global IDs
    DAGCell* c = new DAGCell();
    c->id_ = model::DAG->id_by_index(3, i+1);
    c->dagmc_ptr_ = model::DAG;
    c->universe_ = dagmc_univ_id; // set to zero for now
    c->fill_ = C_NONE; // no fill, single universe

    model::cells.push_back(c);
    model::cell_map[c->id_] = i;

    // Populate the Universe vector and dict
    auto it = model::universe_map.find(dagmc_univ_id);
    if (it == model::universe_map.end()) {
      model::universes.push_back(new Universe());
      model::universes.back()-> id_ = dagmc_univ_id;
      model::universes.back()->cells_.push_back(i);
      model::universe_map[dagmc_univ_id] = model::universes.size() - 1;
    } else {
      model::universes[it->second]->cells_.push_back(i);
    }

    if (model::DAG->is_implicit_complement(vol_handle)) {
      // assuming implicit complement is void for now
      c->material_.push_back(MATERIAL_VOID);
      continue;
    }

    std::string mat_value;
    if (model::DAG->has_prop(vol_handle, "mat")) {
      rval = model::DAG->prop_value(vol_handle, "mat", mat_value);
      MB_CHK_ERR_CONT(rval);
    } else {
      std::stringstream err_msg;
      err_msg << "Volume " << c->id_ << " has no material assignment.";
      fatal_error(err_msg.str());
    }

    std::string cmp_str = mat_value;
    to_lower(cmp_str);
    if (cmp_str.find("void") != std::string::npos   ||
        cmp_str.find("vacuum") != std::string::npos ||
        cmp_str.find("graveyard") != std::string::npos) {
      c->material_.push_back(MATERIAL_VOID);
    } else {
      if (using_uwuw) {
        std::string uwuw_mat = DMD.volume_material_property_data_eh[vol_handle];
        if (uwuw.material_library.count(uwuw_mat) != 0) {
          int matnumber = uwuw.material_library[uwuw_mat].metadata["mat_number"].asInt();
          c->material_.push_back(matnumber);
        } else {
          std::stringstream err_msg;
          err_msg << "Material with value " << mat_value << " not found ";
          err_msg << "in the material library";
          fatal_error(err_msg.str());
        }
      } else {
        c->material_.push_back(std::stoi(mat_value));
      }
    }
  }

  // Allocate the cell overlap count if necessary.
  if (settings::check_overlaps) {
    model::overlap_check_count.resize(model::cells.size(), 0);
  }

  // initialize surface objects
  int n_surfaces = model::DAG->num_entities(2);
  model::surfaces.resize(n_surfaces);

  for (int i = 0; i < n_surfaces; i++) {
    moab::EntityHandle surf_handle = model::DAG->entity_by_index(2, i+1);

    // set cell ids using global IDs
    DAGSurface* s = new DAGSurface();
    s->id_ = model::DAG->id_by_index(2, i+1);
    s->dagmc_ptr_ = model::DAG;

    std::string bc_value;
    if (model::DAG->has_prop(surf_handle, "boundary")) {
      rval = model::DAG->prop_value(surf_handle, "boundary", bc_value);
      MB_CHK_ERR_CONT(rval);
      to_lower(bc_value);    
    
      if (bc_value == "transmit" || bc_value == "transmission") {
        s->bc_ = BC_TRANSMIT;
      } else if (bc_value == "vacuum") {
        s->bc_ = BC_VACUUM;
      } else if (bc_value == "reflective" || bc_value == "reflect" || bc_value == "reflecting") {
        s->bc_ = BC_REFLECT;
      } else if (bc_value == "periodic") {
        fatal_error("Periodic boundary condition not supported in DAGMC.");
      } else {
        std::stringstream err_msg;
        err_msg << "Unknown boundary condition \"" << s->bc_
                << "\" specified on surface " << s->id_;
        fatal_error(err_msg);
      }
    } else { // if no condition is found, set to transmit
      s->bc_ = BC_TRANSMIT;
    }

    // add to global array and map
    model::surfaces[i] = s;
    model::surface_map[s->id_] = s->id_;
  }

  return;
}

void free_memory_dagmc()
{
  delete model::DAG;
}

}
#endif
