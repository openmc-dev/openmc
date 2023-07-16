#ifndef OPENMC_PARTITIONER_UTILS_H
#define OPENMC_PARTITIONER_UTILS_H

#include "aabb.h"
#include "position.h"
#include "random_dist.h"
#include "timer.h"
#include "universe.h"

#include <float.h>

#include <omp.h>
#include <set>
#include <vector>

namespace openmc {

// This works by allocating large pools of nodes
// When we have to allocate nodes during construction, we allocate one from
// these pools instead of using new[]
// NodeT is our type of node
// NUM_CHILDREN_PER_PARENT is the number of children per parent (e.g. octree
// = 8. kd tree = 2) POOL_SIZE controls how many parents can have children
// allocated from here before we have to reallocate
template<class NodeT, size_t NUM_CHILDREN_PER_PARENT, size_t POOL_SIZE = 16384>
class NodeAllocator {
public:
  NodeAllocator() : last_pool_next_index_(0) { omp_init_lock(&lock_); }

  ~NodeAllocator() { omp_destroy_lock(&lock_); }

  NodeT* allocate()
  {
    omp_set_lock(&lock_);

    // If we have no pools or the pool we are currently allocating from is full,
    // let's allocate a new pool
    if (last_pool_next_index_ == POOL_SIZE || pools_.size() == 0) {
      pools_.push_back(
        std::make_unique<NodeT[]>(NUM_CHILDREN_PER_PARENT * POOL_SIZE));

      // Reset the position we allocate nodes from
      last_pool_next_index_ = 0;
    }

    // Get the next abalaible nodes
    auto ptr =
      &pools_.back()[NUM_CHILDREN_PER_PARENT * last_pool_next_index_++];

    omp_unset_lock(&lock_);

    return ptr;
  }

private:
  // Vector of memory pools. The allocator allocates nodes from the last pool
  // until it is full, at which point it allocates a new pool
  std::vector<std::unique_ptr<NodeT[]>> pools_;
  // The position in the last pool were we allocate nodes from
  size_t last_pool_next_index_;
  // Lock to allow safe multithreading
  omp_lock_t lock_;
};

// When building an octree or kd-tree, we first need to sample random points
// within the universe to determine where the universe's cells are. This struct
// is used for that task - it contains the position we sample inside the
// universe and cell we found
struct CellPointUncompressed {
  // The position. Might be better to use single precision for this for better
  // memory efficiency
  Position pos_;
  // The cell ID
  int cell_;

  // Compare cell IDs
  bool operator<(const CellPointUncompressed& other) const;

  // Implicit conversion to int
  operator int() const;
};

// This is an octree-specific struct.
// The uncompressed node takes up a lot of memory for a position variable.
// We often need to quickly iterate through arrays of millions of points, so
// this is bad For octrees. We can observe we only use the position to determine
// which child node the point goes into. We can actually determine which child
// node indicies a point will go into for each depth ahead of time. We can use
// this to our advantage and only store those indicies as 3 bit values in an 64
// bit integer I allocate the least signficant 45 bits to this information (15
// depths) The remaining 19 bits are dedicated to the cell ID
struct CellPoint {
  // The variable we encode our data in
  uint64_t data_;

  // Extract the cell from data_
  int get_cell() const;
  // Extract the child index from data_ for a particular depth
  int get_child_index(int depth) const;

  // Store an uncompressed node into data_
  void compress_from(const CellPointUncompressed& uncomp, const AABB& bounds);

  // Compare cell IDs
  bool operator<(const CellPoint& other) const;
  // Implicit conversion to int
  operator int() const;
};

// During our intital point search, we don't want to uniformly sample points
// because the same number of points will be dedicated to regions with fine
// geometry as with coarse geometry Using a binned point search, we can
// importance sample regions of more complex geometry by learning over time in
// what regions we find a greater variety of cells
class Bin {
public:
  Bin();

  // Insert a cell into the bin (thread safe)
  void insert(int cell);
  // Insert a cell into the bin (not thread safe)
  void insert_lockless(int cell);
  // Sort the cells while neatly encapsulating the Bin class
  void sort_cells();
  // Run an in-place algorithm that makes all the cells in cells_ unique
  void make_cells_unique();
  // Compute the score for this bin, it allows us to assign a probability to it
  // when sampling bins
  double compute_score(int iteration);
  // Get the cells in this bin
  const std::vector<int>& get_cells() const;

  // If we have already sent the bin's cells to Universe::find_cell but have
  // found nothing, we'll have to use the fallback. This allows us to skip the
  // cells we already tested that are in the fallback result
  void copy_untested_cells(
    const std::vector<int>& possible_cells, std::vector<int>& untested_cells);

private:
  // utility functions
  void lock_bin();
  void unlock_bin();
  // Thread lock
  omp_lock_t lock_;
  // Cells in this bin
  std::vector<int> cells_;
  // Number of sorted cells, used in checking whether a bin already has a cell
  int num_sorted_cells_;
  // Variable used in scoring
  int prev_cell_count_;
};

// A bin that's more oriented towards random refinement
template<class NodeT>
class ProbabilityBin {
public:
  ProbabilityBin() : num_searched_(0), num_found_(0)
  {
    omp_init_lock(&found_lock_);
  }

  ~ProbabilityBin()
  {
    omp_destroy_lock(&found_lock_);

    for (auto& p : contained_nodes_) {
      omp_destroy_lock(&p.second);
    }
  }

  // like the regular bin, we have node sampling
  double compute_score()
  {
    // Current, the score is how many nodes are within this bin
    return contained_nodes_.size();
  }

  // get a random node index
  size_t pick_node(uint64_t* seed)
  {
    return (size_t)uniform_distribution(
      0.0, (double)contained_nodes_.size(), seed);
  }

  // Take the unique cells stored in common_cells_ and found_cells_ and store
  // them in the vector
  void update_common_cells()
  {
    // merge the two vectors
    // here I repurpose found_cells_ to avoid allocating memory
    found_cells_.reserve(common_cells_.size() + found_cells_.size());
    for (int cell : common_cells_) {
      found_cells_.push_back(cell);
    }

    // copy over the unique cells from the merged vector into common_cells_
    std::sort(found_cells_.begin(), found_cells_.end());
    common_cells_.clear();
    int prev_cell = -1;
    for (int cell : found_cells_) {
      if (cell != prev_cell) {
        common_cells_.push_back(cell);
        prev_cell = cell;
      }
    }

    // clear the merge vector to allow it to be used as found_cells_ again next
    // iteration
    found_cells_.clear();
  }

  void add_cell(int cell)
  {
    omp_set_lock(&found_lock_);
    found_cells_.push_back(cell);
    omp_unset_lock(&found_lock_);
  }

  int num_searched_, num_found_;
  omp_lock_t found_lock_;
  std::vector<int> common_cells_;
  std::vector<int> found_cells_;
  std::vector<std::pair<NodeT*, omp_lock_t>> contained_nodes_;
};

// Function to automatically store the point into the type given by the vector
// Allows us to use the same binned point search code for octrees and kd trees
// even though octree's binned point search returns a vector of compressed cell
// points whereas kd trees need uncompressed cell points
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

// Return a large point cloud of cells in the scene
template<typename T>
std::vector<T> binned_point_search(
  const Universe& univ, const UniversePartitioner& fallback, const AABB& bounds)
{
  constexpr int32_t BINNING_SEARCH_TOTAL_POINTS =
    8000000; // how many total points we want to search
  constexpr int32_t BINNING_SEARCH_TOTAL_ITERATIONS =
    128; // how many iterations we want to have (many iterations is slightly
         // better but might lead to bad CPU usage)
  constexpr int32_t BINNING_SEARCH_GRID_RES =
    32; // the resolution of our bin grid

  // the density we want to sample per iteration
  double target_density = BINNING_SEARCH_TOTAL_POINTS /
                          (BINNING_SEARCH_TOTAL_ITERATIONS * bounds.volume());

  // the actual bins
  std::vector<Bin> bin_grid(BINNING_SEARCH_GRID_RES * BINNING_SEARCH_GRID_RES *
                            BINNING_SEARCH_GRID_RES);
  // a cumulative mass function that is utilized when sampling bins
  std::vector<double> bin_cmf(BINNING_SEARCH_GRID_RES *
                              BINNING_SEARCH_GRID_RES *
                              BINNING_SEARCH_GRID_RES);

  // size of each bin on each dimension
  Position bin_dim;
  for (int i = 0; i < 3; i++) {
    bin_dim[i] = (bounds.max_[i] - bounds.min_[i]) / BINNING_SEARCH_GRID_RES;
  }

  // account for floating point errors
  int total_search_points =
    BINNING_SEARCH_TOTAL_ITERATIONS * (int)(bounds.volume() * target_density);
  std::vector<T> points_in_bounds(total_search_points);

  // an offset into points_in_bounds marking where we begin to write out the
  // cell points
  int write_offset = 0;
  for (int iteration = 0; iteration < BINNING_SEARCH_TOTAL_ITERATIONS;
       iteration++) {
    int num_search_points = static_cast<int>(bounds.volume() * target_density);

    // init the cmf for the bins
    double total_cmf = 0.0;
    for (int i = 0; i < bin_grid.size(); i++) {
      total_cmf += bin_grid[i].compute_score(iteration);
      bin_cmf[i] = total_cmf;
    }

    for (double& cmf : bin_cmf) {
      cmf /= total_cmf;
    }

#pragma omp parallel
    {
      // initialize per-thread resources
      int tid = omp_get_thread_num();
      int tcount = omp_get_num_threads();

      uint64_t bin_seed = write_offset + tid;
      uint64_t pos_seed[3] = {
        static_cast<uint64_t>(write_offset + tid + tcount),
        static_cast<uint64_t>(write_offset + tid + tcount * 2),
        static_cast<uint64_t>(write_offset + tid + tcount * 3)};

      // temporary buffer
      std::vector<int> untested_cells;
      untested_cells.reserve(univ.cells_.size());

#pragma omp for
      for (int i = 0; i < num_search_points; i++) {
        // index that is written to
        int index = write_offset + i;

        // pick a bin
        double bin_cmf_val = uniform_distribution(0.0, 1.0, &bin_seed);
        auto cmf_iter =
          std::upper_bound(bin_cmf.begin(), bin_cmf.end(), bin_cmf_val);
        int bin_idx = std::distance(bin_cmf.begin(), cmf_iter);
        Bin* sample_bin = &bin_grid[bin_idx];

        // generate a point inside that bin
        CellPointUncompressed point;
        for (int j = 0; j < 3; j++) {
          int idx = bin_idx % BINNING_SEARCH_GRID_RES;
          bin_idx /= BINNING_SEARCH_GRID_RES;
          double min = idx * bin_dim[j] + bounds.min_[j];
          double max = bin_dim[j] + min;

          point.pos_[j] = uniform_distribution(min, max, &pos_seed[j]);
        }

        // check if the bin already has the cell containing our point
        Direction dummy_dir {1.0, 0.0, 0.0};
        point.cell_ =
          univ.find_cell_for_point(sample_bin->get_cells(), point.pos_);

        // if we didn't find the point, we need to sample and reinsert
        if (point.cell_ == -1) {
          const auto& possible_cells =
            fallback.get_cells(point.pos_, dummy_dir);
          sample_bin->copy_untested_cells(possible_cells, untested_cells);
          point.cell_ = univ.find_cell_for_point(untested_cells, point.pos_);

          sample_bin->insert(point.cell_);
        } // else don't bother inserting

        // now store it in the large list of cells that we return to the
        // partitioner
        store_cell_point(points_in_bounds[index], point, bounds);
      }
    }

    // sort the cells afterward
#pragma omp parallel for
    for (auto& bin : bin_grid) {
      bin.sort_cells();
    }

    // move the write offset
    write_offset += num_search_points;
  }

  return points_in_bounds;
}

}; // namespace openmc

#endif // OPENMC_PARTITIONER_UTILS_H