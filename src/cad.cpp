
#ifdef CAD

#include "cad.h"
#include "error.h"

#include <string>
#include <algorithm>

moab::DagMC* DAGMC;

#define TOLOWER(S) std::transform(S.begin(), S.end(), S.begin(), ::tolower)

void load_cad_geometry_c()
{
  if(!DAGMC) {
    DAGMC = new moab::DagMC();
  }

  int32_t cad_univ_id = 0; // universe is always 0 for CAD
  
  moab::ErrorCode rval = DAGMC->load_file("dagmc.h5m");
  MB_CHK_ERR_CONT(rval);

  rval = DAGMC->init_OBBTree();
  MB_CHK_ERR_CONT(rval);

  std::vector< std::string > prop_keywords;
  prop_keywords.push_back("mat");
  prop_keywords.push_back("boundary");
  
  std::map<std::string, std::string> ph;
  DAGMC->parse_properties(prop_keywords, ph, ":");
  MB_CHK_ERR_CONT(rval);
  
  // initialize cell objects
  openmc::n_cells = DAGMC->num_entities(3);
  for(int i = 0; i < openmc::n_cells; i++)
    {
      moab::EntityHandle vol_handle = DAGMC->entity_by_index(3, i+1);
      
      // set cell ids using global IDs
      openmc::CADCell* c = new openmc::CADCell();
      c->id_ = DAGMC->id_by_index(3, i+1);
      c->dagmc_ptr = DAGMC;
      c->universe_ = cad_univ_id; // set to zero for now

      if(DAGMC->has_prop(vol_handle, "mat")){
	std::string mat_value;
	
	rval = DAGMC->prop_value(vol_handle, "mat", mat_value);
	MB_CHK_ERR_CONT(rval);	
	int mat_id = std::stoi(mat_value);
	c->material_.push_back(mat_id);
      }
      else {
	c->material_.push_back(40); // TO-DO: add void material here
      }
      
      c->fill_ = openmc::C_NONE;

      openmc::global_cells.push_back(c);
      openmc::cell_map[c->id_] = c->id_;

      // Populate the Universe vector and dict
      auto it = openmc::universe_map.find(cad_univ_id);
      if (it == openmc::universe_map.end()) {
	openmc::global_universes.push_back(new openmc::Universe());
	openmc::global_universes.back()-> id = cad_univ_id;
	openmc::global_universes.back()->cells.push_back(i);
	openmc::universe_map[cad_univ_id] = openmc::global_universes.size() - 1;
      }
      else {
	openmc::global_universes[it->second]->cells.push_back(i);
      }

      
      if(DAGMC->is_implicit_complement(vol_handle)) {
	// assuming implicit complement is void for now
        c->material.push_back(openmc::MATERIAL_VOID);	
	continue;
      }
      
      if(DAGMC->has_prop(vol_handle, "mat")){
	std::string mat_value;
	rval = DAGMC->prop_value(vol_handle, "mat", mat_value);
	MB_CHK_ERR_CONT(rval);
	TOLOWER(mat_value);
	
	if(mat_value == "void") {
	  c->material.push_back(openmc::MATERIAL_VOID);
	}
	else {
	  c->material.push_back(std::stoi(mat_value));
	}
      }
      else {
	std::cout << "Warning: volume without material found!" << std::endl;
	c->material.push_back(openmc::MATERIAL_VOID);
      }     
    }

  // initialize surface objects
  openmc::n_surfaces = DAGMC->num_entities(2);
  openmc::global_surfaces.resize(openmc::n_surfaces);
  
  for(int i = 0; i < openmc::n_surfaces; i++)
    {
      moab::EntityHandle surf_handle = DAGMC->entity_by_index(2, i+1);
      
      // set cell ids using global IDs
      openmc::CADSurface* s = new openmc::CADSurface();
      s->id = DAGMC->id_by_index(2, i+1);
      s->dagmc_ptr = DAGMC;

      if(DAGMC->has_prop(surf_handle, "boundary")) {
	std::string bc_value;
	rval = DAGMC->prop_value(surf_handle, "boundary", bc_value);
	MB_CHK_ERR_CONT(rval);
	TOLOWER(bc_value);

	if(bc_value == "transmit" || bc_value == "transmission") {
	  s->bc = openmc::BC_TRANSMIT;
	}
	else if(bc_value == "vacuum") {
	  s->bc = openmc::BC_VACUUM;
	}
	else if(bc_value == "reflective" || bc_value == "reflect" || bc_value == "reflecting") {
	  s->bc = openmc::BC_REFLECT;
	}
	else if(bc_value == "periodic") {
	  openmc::fatal_error("Periodic boundary condition not supported in CAD.");
	}
	else {
	  std::stringstream err_msg;
	  err_msg << "Unknown boundary condition \"" << s->bc
		  << "\" specified on surface " << s->id;
	  openmc::fatal_error(err_msg);
	}
      }
      // if no BC property is found, set to transmit
      else {
	s->bc = openmc::BC_TRANSMIT;
      }

      // add to global array and map
      openmc::global_surfaces[i] = s;
      openmc::surface_map[s->id] = s->id;
    }

  return;
}

void free_memory_cad_c()
{
  delete DAGMC;
}

#endif
