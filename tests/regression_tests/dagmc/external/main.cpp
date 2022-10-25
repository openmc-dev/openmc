#include "openmc/capi.h"
#include "openmc/cross_sections.h"
#include "openmc/dagmc.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/geometry_aux.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"
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

  // Create DAGMC ptr
  std::string filename = "dagmc.h5m";
  std::shared_ptr<moab::DagMC> dag_ptr = std::make_shared<moab::DagMC>();
  moab::ErrorCode rval = dag_ptr->load_file(filename.c_str());
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to load file");
  }

  // Initialize acceleration data structures
  rval = dag_ptr->init_OBBTree();
  if (rval != moab::MB_SUCCESS) {
    fatal_error("Failed to initialise OBB tree");
  }

  // Get rid of existing geometry
  std::unordered_map<std::string, int> nuclide_map_copy =
    openmc::data::nuclide_map;
  openmc::data::nuclides.clear();
  openmc::data::nuclide_map = nuclide_map_copy;
  openmc::model::surfaces.clear();
  openmc::model::surface_map.clear();
  openmc::model::cells.clear();
  openmc::model::cell_map.clear();
  openmc::model::universes.clear();
  openmc::model::universe_map.clear();

  // Update materials (emulate an external program)
  for (auto& mat_ptr : openmc::model::materials) {
    mat_ptr->set_temperature(300);
  }

  // Create new DAGMC universe
  openmc::model::universes.push_back(
    std::make_unique<openmc::DAGUniverse>(dag_ptr));
  model::universe_map[model::universes.back()->id_] =
    model::universes.size() - 1;

  // Add cells to universes
  openmc::populate_universes();

  // Set root universe
  openmc::model::root_universe = openmc::find_root_universe();
  openmc::check_dagmc_root_univ();

  // Final geometry setup and assign temperatures
  openmc::finalize_geometry();

  // Finalize cross sections having assigned temperatures
  openmc::finalize_cross_sections();

  // Check that we correctly assigned cell temperatures with non-void fill
  for (auto& cell_ptr : openmc::model::cells) {
    if (cell_ptr->material_.front() != openmc::C_NONE &&
        cell_ptr->temperature() != 300) {
      fatal_error("Failed to set cell temperature");
    }
  }

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
