
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
      c->id_ = DAGMC->id_by_index(3, i);
      c->dagmc_ptr = DAGMC;
      c->universe_ = cad_univ_id; // set to zero for now
      c->material_.push_back(40); // TEMPORARY
      openmc::cells_c.push_back(c);
      openmc::cell_dict[c->id_] = i;

      // Populate the Universe vector and dict
      auto it = openmc::universe_dict.find(cad_univ_id);
      if (it == openmc::universe_dict.end()) {
	openmc::universes_c.push_back(new openmc::Universe());
	openmc::universes_c.back()-> id = cad_univ_id;
	openmc::universes_c.back()->cells.push_back(i);
	openmc::universe_dict[cad_univ_id] = openmc::universes_c.size() - 1;
      }
      else {
	openmc::universes_c[it->second]->cells.push_back(i);
      }
    }

  // initialize surface objects
  openmc::n_surfaces = DAGMC->num_entities(2);
  openmc::surfaces_c = new openmc::Surface*[openmc::n_surfaces];
  
  for(int i = 0; i < openmc::n_surfaces; i++)
    {
      // set cell ids using global IDs
      openmc::CADSurface* s = new openmc::CADSurface();
      s->id = DAGMC->id_by_index(2, i);
      s->dagmc_ptr = DAGMC;
      s->bc = openmc::BC_TRANSMIT;
      openmc::surfaces_c[i] = s;
      openmc::surface_dict[s->id] = i;
    }

  return;
}

void free_memory_cad_c()
{
  delete DAGMC;
}
