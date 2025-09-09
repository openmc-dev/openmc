#ifdef OPENMC_MPI
#include <mpi.h>
#endif

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/geometry_aux.h"
#include "openmc/message_passing.h"
#include "openmc/summary.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_cell.h"
#include "openmc/tallies/tally.h"

using namespace openmc;

int main(int argc, char** argv)
{
#ifdef OPENMC_MPI
  MPI_Comm world {MPI_COMM_WORLD};
  int err = openmc_init(argc, argv, &world);
#else
  int err = openmc_init(argc, argv, nullptr);
#endif
  if (err)
    fatal_error(openmc_err_msg);

  // create a new cell filter
  auto cell_filter = Filter::create<CellFilter>();

  // add all cells to the cell filter
  std::vector<int32_t> cell_indices;
  for (auto& entry : openmc::model::cell_map) {
    cell_indices.push_back(entry.second);
  }
  // enable distribcells offsets for all cells
  prepare_distribcell(&cell_indices);
  // sort to make sure the cell bins appear in the same
  // order as the test relying on the openmc exe
  std::sort(cell_indices.begin(), cell_indices.end());
  cell_filter->set_cells(cell_indices);

  // create a new tally
  auto tally = Tally::create();
  std::vector<Filter*> filters = {cell_filter};
  tally->set_filters(filters);
  tally->set_scores({"flux"});

  // set the temperature of the cell containing
  // the lattice
  auto& root_univ = openmc::model::universes[openmc::model::root_universe];
  auto& lattice_cell = openmc::model::cells[root_univ->cells_[0]];
  lattice_cell->set_temperature(300.0, 0, true);

  // check that material-filled cells return no contained cells
  for (auto& cell : openmc::model::cells) {
    if (cell->type_ == Fill::MATERIAL) {
      auto contained_cells = cell->get_contained_cells();
      assert(contained_cells.empty());
    }
  }

  // set a higher temperature for only one of the lattice cells (ID is 4 in the
  // model)
  model::cells[model::cell_map[4]]->set_temperature(400.0, 3, true);

  // set the density of another lattice cell to 2
  model::cells[model::cell_map[4]]->set_density(2.0, 2, true);

  // the summary file will be used to check that
  // temperatures were set correctly so clear
  // error output can be provided
#ifdef OPENMC_MPI
  if (openmc::mpi::master)
    openmc::write_summary();
#else
  openmc::write_summary();
#endif

  openmc_run();
  openmc_finalize();

#ifdef OPENMC_MPI
  MPI_Finalize();
#endif

  return 0;
}
