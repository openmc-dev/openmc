#include "moab/Core.hpp"
#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/mesh.h"
#include "openmc/tallies/filter_mesh.h"
#include "openmc/tallies/tally.h"
#include <iostream>

int main(int argc, char* argv[])
{

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
  } else {
    std::cout << "Loaded external MOAB mesh from file " << filename
              << std::endl;
  }

  // Add a new unstructured mesh to openmc using new constructor
  model::meshes.push_back(std::make_unique<MOABMesh>(moabPtrLocal));

  // Check we now have 2 copies of the shared ptr
  if (moabPtrLocal.use_count() != 2) {
    fatal_error("Incorrect number of MOAB shared pointers");
  }

  // Auto-assign mesh ID
  model::meshes.back()->set_id(C_NONE);
  int mesh_id = model::meshes.back()->id_;

  // Check we now have 2 meshes and id was correctly set
  if (model::meshes.size() != 2)
    fatal_error("Wrong number of meshes.");
  else if (mesh_id != 2)
    fatal_error("Mesh ID is incorrect");

  // Add a new mesh filter with auto-assigned ID
  Filter* filter_ptr = Filter::create("mesh", C_NONE);

  // Upcast pointer type
  MeshFilter* mesh_filter = dynamic_cast<MeshFilter*>(filter_ptr);

  if (!mesh_filter) {
    fatal_error("Failed to create mesh filter");
  }

  // Pass in the index of our mesh to the filter
  int32_t mesh_idx = model::meshes.size() - 1;
  mesh_filter->set_mesh(mesh_idx);

  // Create a tally with auto-assigned ID
  model::tallies.push_back(make_unique<Tally>(C_NONE));

  // Set tally name - matches that in test.py
  model::tallies.back()->name_ = "external mesh tally";

  // Set tally filter to our mesh filter
  std::vector<Filter*> filters(1, filter_ptr);
  model::tallies.back()->set_filters(filters);

  // Set tally estimator
  model::tallies.back()->estimator_ = TallyEstimator::TRACKLENGTH;

  // Set tally score
  std::vector<std::string> score_names(1, "flux");
  model::tallies.back()->set_scores(score_names);

  // Run OpenMC
  openmc_err = openmc_run();
  if (openmc_err)
    fatal_error(openmc_err_msg);

  // Deallocate memory
  openmc_err = openmc_finalize();
  if (openmc_err)
    fatal_error(openmc_err_msg);

  return EXIT_SUCCESS;
}
