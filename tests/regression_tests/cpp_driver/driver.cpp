
#include "openmc/capi.h"
#include "openmc/cell.h"
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

  openmc_run();
  openmc_finalize();
  return 0;
}
