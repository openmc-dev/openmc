#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/mesh.h"
#include <iostream>
#include "moab/Core.hpp"

int main(int argc, char * argv[]){

  using namespace openmc;
  int openmc_err;

  // Initialise OpenMC
  openmc_err = openmc_init(argc, argv, nullptr);
  if (openmc_err == -1) {
    // This happens for the -h and -v flags
    return EXIT_SUCCESS;
  } else if (openmc_err) {
    fatal_error(openmc_err_msg);
  }

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
    std::cout<<"Loaded external MOAB mesh from file "<<filename<<std::endl;
  }

  // Add a new unstructured mesh to openmc using new constructor
  model::meshes.push_back(std::make_unique<MOABMesh>(moabPtrLocal));

  // Check we now have 2 copies of the shared ptr
  if(moabPtrLocal.use_count()!=2){
     fatal_error("Incorrect number of MOAB shared pointers");
  }

  // Set mesh ID (-1 signifies auto-assign)
  model::meshes.back()->set_id(-1);
  int mesh_id = model::meshes.back()->id_;

  // Check we now have 2 meshes and id was correctly set
  if(model::meshes.size()!=2)
     fatal_error("Wrong number of meshes.");
   else if(mesh_id!=2)
     fatal_error("Mesh ID is incorrect");

  // Add a new mesh filter tally

  // Run OpenMC
  openmc_err = openmc_run();
  if (openmc_err) fatal_error(openmc_err_msg);

  // Deallocate memory
  openmc_err = openmc_finalize();
    if (openmc_err) fatal_error(openmc_err_msg);

  return EXIT_SUCCESS;

}
