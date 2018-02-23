
#include "cad.h"

extern moab::DagMC* DAGMC;

namespace moab {

void load_cad_geometry_c()
{
  if(!DAGMC) {
    DAGMC = new DagMC();
  }

  ErrorCode rval = DAGMC->load_file("dagmc.h5m");
  MB_CHK_ERR_CONT(rval);

  rval = DAGMC->init_OBBTree();
  MB_CHK_ERR_CONT(rval);

  // initialize cell objects
  for(int i = 0; i < DAGMC->num_entities(3); i++){

    // set cell ids using global IDs
    openmc::CADCell c = openmc::CADCell();
    c.id = DAGMC->id_by_index(3, 0);
    openmc::cells_c.push_back(c);
  }
  
  return;
}
  
}
