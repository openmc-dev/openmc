
#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/geometry.h"
#include "openmc/summary.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_cell.h"
#include "openmc/tallies/tally.h"

using namespace openmc;

int main(int argc, char** argv) {
  openmc_init(argc, argv, nullptr);

  // create a new cell filter
  auto cell_filter = Filter::create<CellFilter>();

  // add all cells to the cell filter
  std::vector<int32_t> cell_indices;
  for (auto& entry : openmc::model::cell_map) {
      cell_indices.push_back(entry.second);
  }
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
  lattice_cell->set_temperature(300, 1, true);

  // check that material-filled cells return no contained cells
  for (auto& cell : openmc::model::cells) {
    if (cell->type_ == Fill::MATERIAL) {
      auto contained_cells = cell->get_contained_cells();
      assert(contained_cells.empty());
    }
  }

  // the summary file will be used to check that
  // temperatures were set correctly so clear
  // error output can be provided
  openmc::write_summary();

  openmc_run();
  openmc_finalize();
  return 0;
}
