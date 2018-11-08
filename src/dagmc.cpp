#include "openmc/dagmc.h"

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/string_utils.h"
#include "openmc/settings.h"
#include "openmc/geometry.h"

#include <string>
#include <sstream>
#include <algorithm>

#ifdef DAGMC

namespace openmc {

moab::DagMC* DAG;

void load_dagmc_geometry()
{
  if (!DAG) {
    DAG = new moab::DagMC();
  }

  int32_t dagmc_univ_id = 0; // universe is always 0 for DAGMC

  moab::ErrorCode rval = DAG->load_file("dagmc.h5m");
  MB_CHK_ERR_CONT(rval);

  rval = DAG->init_OBBTree();
  MB_CHK_ERR_CONT(rval);

  std::vector<std::string> prop_keywords;
  prop_keywords.push_back("mat");
  prop_keywords.push_back("boundary");

  std::map<std::string, std::string> ph;
  DAG->parse_properties(prop_keywords, ph, ":");
  MB_CHK_ERR_CONT(rval);

  // initialize cell objects
  model::n_cells = DAG->num_entities(3);

  // Allocate the cell overlap count if necessary.
  if (settings::check_overlaps) overlap_check_count.resize(model::n_cells, 0);

  for (int i = 0; i < model::n_cells; i++) {
    moab::EntityHandle vol_handle = DAG->entity_by_index(3, i+1);

    // set cell ids using global IDs
    DAGCell* c = new DAGCell();
    c->id_ = DAG->id_by_index(3, i+1);
    c->dagmc_ptr_ = DAG;
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

    if (DAG->is_implicit_complement(vol_handle)) {
      // assuming implicit complement is void for now
      c->material_.push_back(MATERIAL_VOID);
      continue;
    }

    if (DAG->has_prop(vol_handle, "mat")){
      std::string mat_value;
      rval = DAG->prop_value(vol_handle, "mat", mat_value);
      MB_CHK_ERR_CONT(rval);
      to_lower(mat_value);

      if (mat_value == "void" || mat_value == "vacuum") {
        c->material_.push_back(MATERIAL_VOID);
      } else {
        c->material_.push_back(std::stoi(mat_value));
      }
    } else {
      std::stringstream err_msg;
      err_msg << "Volume " << c->id_ << " has no material assignment.";
      fatal_error(err_msg.str());
    }
  }

  // initialize surface objects
  n_surfaces = DAG->num_entities(2);
  surfaces.resize(n_surfaces);

  for (int i = 0; i < n_surfaces; i++) {
    moab::EntityHandle surf_handle = DAG->entity_by_index(2, i+1);

    // set cell ids using global IDs
    DAGSurface* s = new DAGSurface();
    s->id_ = DAG->id_by_index(2, i+1);
    s->dagmc_ptr_ = DAG;

    if (DAG->has_prop(surf_handle, "boundary")) {
      std::string bc_value;
      rval = DAG->prop_value(surf_handle, "boundary", bc_value);
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
    } else {   // if no BC property is found, set to transmit
      s->bc_ = BC_TRANSMIT;
    }

    // add to global array and map
    surfaces[i] = s;
    surface_map[s->id_] = s->id_;
  }

  return;
}

void free_memory_dagmc()
{
  delete DAG;
}

}
#endif
