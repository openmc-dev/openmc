
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

  // initialize cell objects
  openmc::n_cells = DAGMC->num_entities(3);
  for(int i = 0; i < openmc::n_cells; i++)
    {
      // set cell ids using global IDs
      openmc::CADCell c = openmc::CADCell();
      c.id = DAGMC->id_by_index(3, i);
      c.dagmc_ptr = DAGMC;
      openmc::cells_c.push_back(c);
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
      openmc::surfaces_c[i] = s;
    }

  return;
}

void free_memory_cad_c()
{
  delete DAGMC;
}
