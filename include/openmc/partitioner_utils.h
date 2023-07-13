#ifndef OPENMC_PARTITIONER_UTILS_H
#define OPENMC_PARTITIONER_UTILS_H

#include "position.h"
#include "random_dist.h"
#include "timer.h"
#include "universe.h"

#include <float.h>

#include <omp.h>
#include <set>
#include <vector>

namespace openmc {

// This axis-aligned bounding box class (AABB) is designed to work easier with
// partitioners than openmc::BoundingBox
struct AABB {
  AABB();
  AABB(const Position& mi, const Position& ma);

  void extend_max(const Position& val);
  void extend_min(const Position& val);

  void extend(const Position& pos);
  void extend(const AABB& other_box);

  Position get_center() const;
  double surface_area() const;
  double volume() const;
  bool contains(const Position& pos) const;

  void reset();

  bool operator==(const AABB& other) const;
  bool operator!=(const AABB& other) const;

  Position min;
  Position max;
};

template<class NodeT, size_t NUM_CHILDREN_PER_PARENT, size_t POOL_SIZE = 16384>
class NodeAllocator {
public:
  NodeAllocator() : last_pool_next_index(0) {}

  ~NodeAllocator()
  {
    for (auto* ptr : pools) {
      delete[] ptr;
    }
  }

  NodeT* allocate()
  {
    if (last_pool_next_index == POOL_SIZE || pools.size() == 0) {
      pools.push_back(new NodeT[NUM_CHILDREN_PER_PARENT * POOL_SIZE]);
      last_pool_next_index = 0;
    }

    return &pools.back()[NUM_CHILDREN_PER_PARENT * last_pool_next_index++];
  }

private:
  std::vector<NodeT*> pools;
  size_t last_pool_next_index;
};

struct CellPointUncompressed {
  Position pos;
  int cell;

  bool operator<(const CellPointUncompressed& other) const;
  operator int() const;
};

struct CellPoint {
  uint64_t data;

  int get_cell() const;
  int get_child_index(int depth) const;
  void compress_from(const CellPointUncompressed& uncomp, const AABB& bounds);

  bool operator<(const CellPoint& other) const;
  operator int() const;
};

class Bin {
public:
  Bin();

  void insert(int cell);
  void insert_lockless(int cell);

  int size() const;
  int unique_size() const;

  void make_cells_unique();
  void sort_cells();

  double score(int iteration);

  const std::vector<int>& get_cells() const;

  void copy_untested_cells(
    const std::vector<int>& possible_cells, std::vector<int>& untested_cells);

private:
  void lock_bin();
  void unlock_bin();

  omp_lock_t lock;
  std::vector<int> cells;
  int num_unique_cells;
  int prev_cell_count;
};

template<class NodeT>
class ProbabilityBin {
public:
  ProbabilityBin() : num_searched(0), num_found(0)
  {
    omp_init_lock(&found_lock);
  }

  ~ProbabilityBin()
  {
    omp_destroy_lock(&found_lock);

    for (auto& p : contained_nodes) {
      omp_destroy_lock(&p.second);
    }
  }

  double compute_score()
  {
    const double REFINEMENT_BIN_SCORE_FOUND_MULT = 1.0;
    const double REFINEMENT_BIN_SCORE_SEARCHED_MULT = 1.0;

    double found_score =
      std::max(REFINEMENT_BIN_SCORE_FOUND_MULT * num_found, 1.0) /
      std::max(REFINEMENT_BIN_SCORE_SEARCHED_MULT * num_searched, 1.0);

    found_score *= found_score;

    double size_score = contained_nodes.size();
    return size_score;
  }

  size_t pick_node(uint64_t* seed)
  {
    return (size_t)uniform_distribution(
      0.0, (double)contained_nodes.size(), seed);
  }

  void update_common_cells()
  {
    found_cells.reserve(common_cells.size() + found_cells.size());

    for (int cell : common_cells) {
      found_cells.push_back(cell);
    }

    std::sort(found_cells.begin(), found_cells.end());

    common_cells.clear();

    int prev_cell = -1;
    for (int cell : found_cells) {
      if (cell != prev_cell) {
        common_cells.push_back(cell);
        prev_cell = cell;
      }
    }

    found_cells.clear();
  }

  void add_cell(int cell)
  {
    omp_set_lock(&found_lock);
    found_cells.push_back(cell);
    omp_unset_lock(&found_lock);
  }

  int num_searched, num_found;
  omp_lock_t found_lock;
  std::vector<int> common_cells;
  std::vector<int> found_cells;
  std::vector<std::pair<NodeT*, omp_lock_t>> contained_nodes;
};

template<typename StoreT, typename ToStoreT>
void store_cell_point(
  StoreT& storage, const ToStoreT& to_store, const AABB& bounds);

template<>
inline void store_cell_point(
  CellPoint& storage, const CellPointUncompressed& to_store, const AABB& bounds)
{
  storage.compress_from(to_store, bounds);
}

template<>
inline void store_cell_point(CellPointUncompressed& storage,
  const CellPointUncompressed& to_store, const AABB& bounds)
{
  storage = to_store;
}

template<typename T>
std::vector<T> binned_point_search(
  const Universe& univ, const UniversePartitioner& fallback, const AABB& bounds)
{
  const int32_t BINNING_SEARCH_TOTAL_POINTS = 8000000;
  const int32_t BINNING_SEARCH_TOTAL_ITERATIONS = 128;
  const int32_t BINNING_SEARCH_GRID_RES = 32;

  double target_density = BINNING_SEARCH_TOTAL_POINTS /
                          (BINNING_SEARCH_TOTAL_ITERATIONS * bounds.volume());

  std::vector<double> search_densities;
  for (int i = 0; i < BINNING_SEARCH_TOTAL_ITERATIONS; i++) {
    search_densities.push_back(target_density);
  }

  // some console output vars
  int total_search_points = 0;
  for (double density : search_densities) {
    total_search_points += (int)(bounds.volume() * density);
  }

  std::vector<Bin> bin_grid(BINNING_SEARCH_GRID_RES * BINNING_SEARCH_GRID_RES *
                            BINNING_SEARCH_GRID_RES);
  std::vector<double> bin_cdf(BINNING_SEARCH_GRID_RES *
                              BINNING_SEARCH_GRID_RES *
                              BINNING_SEARCH_GRID_RES);

  Position bin_dim;
  for (int i = 0; i < 3; i++) {
    bin_dim[i] = (bounds.max[i] - bounds.min[i]) / BINNING_SEARCH_GRID_RES;
  }

  std::vector<T> points_in_bounds(total_search_points);
  int write_offset = 0;
  for (int iteration = 0; iteration < search_densities.size(); iteration++) {
    double density = search_densities[iteration];
    int num_search_points = static_cast<int>(bounds.volume() * density);

    double total_cdf = 0.0;
    for (int i = 0; i < bin_grid.size(); i++) {
      total_cdf += bin_grid[i].score(iteration);
      bin_cdf[i] = total_cdf;
    }

    for (double& cdf : bin_cdf) {
      cdf /= total_cdf;
    }

#pragma omp parallel
    {
      int tid = omp_get_thread_num();
      int tcount = omp_get_num_threads();

      uint64_t bin_seed = write_offset + tid;
      uint64_t pos_seed[3] = {
        static_cast<uint64_t>(write_offset + tid + tcount),
        static_cast<uint64_t>(write_offset + tid + tcount * 2),
        static_cast<uint64_t>(write_offset + tid + tcount * 3)};

      std::vector<int> untested_cells;
      untested_cells.reserve(univ.cells_.size());

#pragma omp for
      for (int i = 0; i < num_search_points; i++) {

        int index = write_offset + i;
        CellPointUncompressed point;

        double bin_cdf_val = uniform_distribution(0.0, 1.0, &bin_seed);
        auto cdf_iter =
          std::upper_bound(bin_cdf.begin(), bin_cdf.end(), bin_cdf_val);
        int bin_idx = std::distance(bin_cdf.begin(), cdf_iter);
        Bin* sample_bin = &bin_grid[bin_idx];

        for (int j = 0; j < 3; j++) {
          int idx = bin_idx % BINNING_SEARCH_GRID_RES;
          bin_idx /= BINNING_SEARCH_GRID_RES;
          double min = idx * bin_dim[j] + bounds.min[j];
          double max = bin_dim[j] + min;

          point.pos[j] = uniform_distribution(min, max, &pos_seed[j]);
        }

        Direction dummy_dir {1.0, 0.0, 0.0};

        point.cell =
          univ.find_cell_for_point(sample_bin->get_cells(), point.pos);
        if (point.cell == -1) {
          const auto& possible_cells = fallback.get_cells(point.pos, dummy_dir);
          sample_bin->copy_untested_cells(possible_cells, untested_cells);
          point.cell = univ.find_cell_for_point(untested_cells, point.pos);

          sample_bin->insert(point.cell);
        } // else don't bother inserting

        store_cell_point(points_in_bounds[index], point, bounds);
      }
    }

#pragma omp parallel for
    for (auto& bin : bin_grid) {
      bin.sort_cells();
    }

    write_offset += num_search_points;
  }

  return points_in_bounds;
}

}; // namespace openmc

#endif // OPENMC_PARTITIONER_UTILS_H