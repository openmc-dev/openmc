#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/mesh.h"
#include <iostream>
#include "moab/Core.hpp"

namespace openmc::model {
  extern std::map<int32_t, std::shared_ptr<moab::Interface> > moabPtrs;
}


int main(int argc, char * argv[]){

  using namespace openmc;
  int openmc_err;

  // Create MOAB interface
  std::shared_ptr<moab::Interface> moabPtrLocal =
    std::make_shared<moab::Core>();

  // Load unstructured mesh file
  std::string filename = "test_mesh_tets.h5m";
  moab::ErrorCode rval = moabPtrLocal->load_file(filename.c_str());
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to load the unstructured mesh file: " + filename);
  }
  else{
    std::cout<<"Loaded external mesh from file "<<filename<<std::endl;
  }

  // Set OpenMC's copy of MOAB (requires fore-knowledge of mesh ID)
  int32_t meshID = 1;
  model::moabPtrs[meshID] = moabPtrLocal;

  // Initialise OpenMC
  openmc_err = openmc_init(argc, argv, nullptr);

  if (openmc_err == -1) {
    // This happens for the -h and -v flags
    return EXIT_SUCCESS;
  } else if (openmc_err) {
    fatal_error(openmc_err_msg);
  }

  // Run OpenMC
  openmc_err = openmc_run();
  if (openmc_err) fatal_error(openmc_err_msg);

  // Deallocate memory
  openmc_err = openmc_finalize();
    if (openmc_err) fatal_error(openmc_err_msg);

  return EXIT_SUCCESS;

}
