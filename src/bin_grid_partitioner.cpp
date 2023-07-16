#include "openmc/bin_grid_partitioner.h"
#include "openmc/error.h"
#include "openmc/partitioner_utils.h"
#include "openmc/timer.h"

#include <set>

namespace openmc {

BinGridPartitioner::BinGridPartitioner(
  const Universe& univ, const AABB& bounds, uint32_t grid_res)
  : grid_res_(grid_res), fallback_(univ), bounds_(bounds)
{
  write_message("Building bin grid partitioner...", 5);
  Timer build_timer;
  build_timer.start();

  // Initialize
  for (int i = 0; i < 3; i++) {
    bin_dim_[i] = (bounds_.max_[i] - bounds_.min_[i]) / grid_res;
  }

  // Put all points into the bins
  auto points =
    binned_point_search<CellPointUncompressed>(univ, fallback_, bounds_);
  bin_grid_.resize(grid_res * grid_res * grid_res);
  for (const auto& p : points) {
    bin_grid_[convert_to_index(p.pos_)].push_back(p.cell_);
  }

  // Remove any redundant points that may be in the bins
  for (auto& bin : bin_grid_) {
    std::set<int> unique_cells;
    for (int cell : bin) {
      unique_cells.insert(cell);
    }
    bin.clear();
    for (int cell : unique_cells) {
      bin.push_back(cell);
    }
  }

  write_message("Bin grid partitioner construction completed in " +
                  std::to_string(build_timer.elapsed()) + " seconds.",
    5);
}

BinGridPartitioner::~BinGridPartitioner() {}

int BinGridPartitioner::convert_to_index(const Position& r) const
{
  // The idea behind this mapping is to count, on each axis, how many bins
  // (rounded down) the point is from the minimum bounds of the bin grid
  // partitioner. This information gives us the location of a unique bins, and
  // we can convert this information into a base grid_res_ number to map it to
  // an index.
  int idx = 0;
  for (int i = 0; i < 3; i++) {
    idx = grid_res_ * idx +
          static_cast<int>((r[i] - bounds_.min_[i]) / bin_dim_[i]);
  }
  return idx;
}

const std::vector<int>& BinGridPartitioner::get_cells(
  Position r, Direction u) const
{
  if (!bounds_.contains(r)) {
    return get_cells_fallback(r, u);
  }

  return bin_grid_[convert_to_index(r)];
}

const std::vector<int>& BinGridPartitioner::get_cells_fallback(
  Position r, Direction u) const
{
  return fallback_.get_cells(r, u);
}

}; // namespace openmc