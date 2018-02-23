
#include "cad.h"

void load_cad_geometry_c()
{
  moab::DagMC* DAGMC = new moab::DagMC();

  moab::ErrorCode rval = DAGMC->load_file("dagmc.h5m");
  MB_CHK_ERR_CONT(rval);

  rval = DAGMC->init_OBBTree();
  MB_CHK_ERR_CONT(rval);

  // initialize cell objects
  for(int i = 0; i < DAGMC->num_entities(3); i++)
    {
      // set cell ids using global IDs
      openmc::CADCell c = openmc::CADCell();
      c.id = DAGMC->id_by_index(3, i);
      openmc::cells_c.push_back(c);
    }
  
  return;
}
