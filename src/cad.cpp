
#include "cad.h"

moab::DagMC* DAGMC;

void load_cad_geometry_c()
{
  if(!DAGMC) {
    DAGMC = new moab::DagMC();
  }

  moab::ErrorCode rval = DAGMC->load_file("dagmc.h5m");
  MB_CHK_ERR_CONT(rval);

  rval = DAGMC->init_OBBTree();
  MB_CHK_ERR_CONT(rval);

  int32_t cad_univ_id = 0; // universe is always 0 for CAD
  
  // initialize cell objects
  openmc::n_cells = DAGMC->num_entities(3);
  for(int i = 0; i < openmc::n_cells; i++)
    {
      // set cell ids using global IDs
      openmc::CADCell* c = new openmc::CADCell();
      c->id_ = DAGMC->id_by_index(3, i+1);
      c->dagmc_ptr = DAGMC;
      c->universe_ = cad_univ_id; // set to zero for now
      c->material_.push_back(40); // TEMPORARY
      openmc::global_cells.push_back(c);
      openmc::cell_map[c->id_] = c->id;

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
    }

  // initialize surface objects
  openmc::n_surfaces = DAGMC->num_entities(2);
  openmc::surfaces_c = new openmc::Surface*[openmc::n_surfaces];
  
  for(int i = 0; i < openmc::n_surfaces; i++)
    {
      // set cell ids using global IDs
      openmc::CADSurface* s = new openmc::CADSurface();
      s->id = DAGMC->id_by_index(2, i+1);
      s->dagmc_ptr = DAGMC;
      s->bc = openmc::BC_TRANSMIT;
      openmc::surfaces_c[i] = s;
      openmc::surface_map[s->id] = s->id;
    }

  return;
}

void free_memory_cad_c()
{
  delete DAGMC;
}
