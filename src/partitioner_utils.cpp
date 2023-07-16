#include "openmc/partitioner_utils.h"
#include <algorithm>
#include <stdint.h>

namespace openmc {

// Compare by cell IDs
bool CellPointUncompressed::operator<(const CellPointUncompressed& other) const
{
  return (cell_ < other.cell_);
}

// Implicit conversion to int
CellPointUncompressed::operator int() const
{
  return cell_;
}

// Store octree child information for only 15 depths (45 bits)
constexpr int32_t CELL_POINT_COMPRESSION_MAX_DEPTH = 15;

// Extract the cell ID
int CellPoint::get_cell() const
{
  uint64_t index = (data_ >> (3 * CELL_POINT_COMPRESSION_MAX_DEPTH));
  return (int)index;
}

// Extract the child index for a particular depth
int CellPoint::get_child_index(int depth) const
{
  uint64_t index = ((data_ >> (3 * depth)) & 0b111);
  return (int)index;
}

// Pack an uncompressed cell
void CellPoint::compress_from(
  const CellPointUncompressed& uncomp, const AABB& bounds)
{
  // Clear data
  data_ = 0;

  // Initialize variables for our encoding loop
  Position node_dim;
  for (int i = 0; i < 3; i++) {
    node_dim[i] = (bounds.max_[i] - bounds.min_[i]) * 0.5;
  }

  // Go through each depth and find what child index we would pick
  // center stores the center of the parent node, starting at the root's center,
  // which is the middle of the bounding box Over time, we update center to find
  // out center of the child node this point would go into
  Position center = bounds.get_center();
  for (int i = 0; i < CELL_POINT_COMPRESSION_MAX_DEPTH; i++) {
    // idx is what child index we would pick
    uint64_t idx = 0;

    for (int i = 0; i < 3; i++) {
      // Determine where the point lies relative to the current parent node
      bool less = (uncomp.pos_[i] < center[i]);

      // Update next center location accordingly
      node_dim[i] *= 0.5;
      center[i] += node_dim[i] * (less ? -1 : 1);

      // Update child index according
      idx = 2 * idx + int(less);
    }

    // Store the index in our data
    data_ |= (idx << (3 * i));
  }

  // Move our cell ID ahead of the bits used to encode the child index
  // information
  data_ |= (static_cast<uint64_t>(uncomp.cell_)
            << (3 * CELL_POINT_COMPRESSION_MAX_DEPTH));
}

// Compare by cell IDs
bool CellPoint::operator<(const CellPoint& other) const
{
  return (get_cell() < other.get_cell());
}

// Implicit conversion to int
CellPoint::operator int() const
{
  return get_cell();
}

// Initialize some variables and reserve memory
Bin::Bin() : num_sorted_cells_(0), prev_cell_count_(0)
{
  omp_init_lock(&lock_);
  cells_.reserve(512);
}

void Bin::insert(int cell)
{
  lock_bin();
  insert_lockless(cell);
  unlock_bin();
}

void Bin::insert_lockless(int cell)
{
  cells_.push_back(cell);
  if (cells_.size() == cells_.capacity()) {
    int cur_size = cells_.size();

    make_cells_unique();

    int size_saving = cur_size - cells_.size();
    if (size_saving < (int)(0.25 * cur_size)) {
      cells_.reserve(2 * cells_.capacity());
    }
  }
}

void Bin::sort_cells()
{
  std::sort(cells_.begin(), cells_.end());
}

void Bin::make_cells_unique()
{
  // First sort the cells
  sort_cells();
  // Compare adjacent cell IDs to find the unique ones and replace the common
  // ones in place
  int i = 0, j = 1;
  while (j < cells_.size()) {
    if (cells_[j] != cells_[i]) {
      i++;
      cells_[i] = cells_[j];
    }
    j++;
  }
  // Make the size smaller so the vector only consists of unique cells
  cells_.resize(i + 1);
}

double Bin::compute_score(int iteration)
{
  // The score we return at the end
  double score;

  // we use max here to prevent cases where we don't sample a bin at all
  size_t bounded_size = std::max(cells_.size(), static_cast<size_t>(1));
  if (iteration == 0) {
    // In the start, we don't know anything about the density of unique cells
    // Let's sample everything uniformly
    score = 1.0;
  } else if (iteration == 1) {
    // We know some information about the scene, so let's set the score to how
    // many unique cells there are
    score = (double)bounded_size;
  } else {
    // Since there have been quite a few iterations, now we can incorporate how
    // many new cells we are finding in a particular bin into the mix If we have
    // exhausted most cells in a particular bin, then we can slighly discourage
    // from resampling it again If we are making many discoveries in another
    // bin, we can put more sampling effort towards it
    constexpr double alpha = 0.1;

    double size_score = bounded_size;
    double cell_increase_score = double(cells_.size()) / prev_cell_count_;

    score = alpha * cell_increase_score + (1 - alpha) * size_score;
  }

  // do stuff here to update for the next scoring cycle
  prev_cell_count_ = bounded_size;

  return score;
}

const std::vector<int>& Bin::get_cells() const
{
  return cells_;
}

// assumes that cells are already sorted
void Bin::copy_untested_cells(
  const std::vector<int>& possible_cells, std::vector<int>& untested_cells)
{
  untested_cells.clear();

  lock_bin();
  // Copy over the cells that are not in the vector
  for (int cell : possible_cells) {
    if (!std::binary_search(
          cells_.begin(), cells_.begin() + num_sorted_cells_, cell)) {
      untested_cells.push_back(cell);
    }
  }
  unlock_bin();
}
void Bin::lock_bin()
{
  omp_set_lock(&lock_);
}

void Bin::unlock_bin()
{
  omp_unset_lock(&lock_);
}

}; // namespace openmc