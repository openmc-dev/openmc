#include "openmc/partitioner_utils.h"
#include <algorithm>
#include <float.h>
#include <stdint.h>

namespace openmc {

const int32_t CELL_POINT_COMPRESSION_MAX_DEPTH = 15;

inline vec3 vmin(vec3 lhs, vec3 rhs)
{
  vec3 res;
  res.x = std::min(lhs.x, rhs.x);
  res.y = std::min(lhs.y, rhs.y);
  res.z = std::min(lhs.z, rhs.z);
  return res;
}

inline vec3 vmax(vec3 lhs, vec3 rhs)
{
  vec3 res;
  res.x = std::max(lhs.x, rhs.x);
  res.y = std::max(lhs.y, rhs.y);
  res.z = std::max(lhs.z, rhs.z);
  return res;
}

AABB::AABB() : min(FLT_MAX, FLT_MAX, FLT_MAX), max(-FLT_MAX, -FLT_MAX, -FLT_MAX)
{}

AABB::AABB(const vec3& mi, const vec3& ma) : min(mi), max(ma) {}

void AABB::reset()
{
  AABB clear;
  *this = clear;
}

void AABB::extend_max(const vec3& val)
{
  max = vmax(max, val);
}

void AABB::extend_min(const vec3& val)
{
  min = vmin(min, val);
}

void AABB::extend(const vec3& pos)
{
  extend_max(pos);
  extend_min(pos);
}

float AABB::surface_area() const
{
  vec3 side_lengths = max - min;

  return 2 * (side_lengths.x * (side_lengths.y + side_lengths.z) +
               side_lengths.y * side_lengths.z);
}

float AABB::volume() const
{
  vec3 side_lengths = max - min;
  return side_lengths.x * side_lengths.y * side_lengths.z;
}

void AABB::extend(const AABB& other_box)
{
  extend_max(other_box.max);
  extend_min(other_box.min);
}

vec3 AABB::get_center() const
{
  return (min + max) * 0.5f;
}

bool AABB::contains(const vec3& pos) const
{
  return (min.x <= pos.x && min.y <= pos.y && min.z <= pos.z &&
          max.x >= pos.x && max.y >= pos.y && max.z >= pos.z);
}

bool AABB::operator==(const AABB& other) const
{
  return (min.x == other.min.x && min.y == other.min.y &&
          min.z == other.min.z && max.x == other.max.x &&
          max.y == other.max.y && max.z == other.max.z);
}

bool AABB::operator!=(const AABB& other) const
{
  return !(*this == other);
}

bool CellPointUncompressed::operator<(const CellPointUncompressed& other) const
{
  return (cell < other.cell);
}

int CellPoint::get_cell() const
{
  uint64_t index = (data >> (3 * CELL_POINT_COMPRESSION_MAX_DEPTH));
  return (int)index;
}

int CellPoint::get_child_index(int depth) const
{
  uint64_t index = ((data >> (3 * depth)) & 0b111);
  return (int)index;
}

void CellPoint::compress_from(
  const CellPointUncompressed& uncomp, const AABB& bounds)
{
  data = 0;

  vec3 node_dim;
  for (int i = 0; i < 3; i++) {
    node_dim[i] = (bounds.max[i] - bounds.min[i]) * 0.5;
  }

  vec3 center = bounds.get_center();

  for (int i = 0; i < CELL_POINT_COMPRESSION_MAX_DEPTH; i++) {
    // halve the node dim
    uint64_t idx = 0;
    for (int i = 0; i < 3; i++) {
      node_dim[i] *= 0.5;
      bool less = (uncomp.pos[i] < center[i]);
      center[i] += node_dim[i] * (less ? -1 : 1);
      idx = 2 * idx + int(less);
    }

    idx = (idx << (3 * i));
    data |= idx;
  }

  uint64_t cell_comp = uncomp.cell;
  cell_comp = (cell_comp << (3 * CELL_POINT_COMPRESSION_MAX_DEPTH));
  data |= cell_comp;
}

bool CellPoint::operator<(const CellPoint& other) const
{
  return (get_cell() < other.get_cell());
}

Bin::Bin() : num_unique_cells(0), prev_cell_count(0)
{
  omp_init_lock(&lock);
  cells.reserve(512);
}

void Bin::insert(int cell)
{
  lock_bin();
  insert_lockless(cell);
  unlock_bin();
}

void Bin::insert_lockless(int cell)
{
  cells.push_back(cell);
  if (cells.size() == cells.capacity()) {
    int cur_size = cells.size();

    make_cells_unique();

    int size_saving = cur_size - cells.size();
    if (size_saving < (int)(0.25 * cur_size)) {
      cells.reserve(2 * cells.capacity());
    }
  }
}

void Bin::sort_cells()
{
  std::sort(cells.begin(), cells.end());
}

void Bin::make_cells_unique()
{
  sort_cells();
  int i = 0, j = 1;
  while (j < cells.size()) {
    if (cells[j] != cells[i]) {
      i++;
      cells[i] = cells[j];
    }
    j++;
  }
  cells.resize(i + 1);
  num_unique_cells = cells.size();
}

int Bin::size() const
{
  return std::max(cells.size(), (size_t)1);
}

int Bin::unique_size() const
{
  return num_unique_cells;
}

float Bin::score(int iteration)
{
  float val;

  if (iteration == 0) {
    val = 1.0;
  } else if (iteration == 1) {
    val = (float)size();
  } else {
    const float alpha = 0.1;

    float size_score = size();
    float cell_increase_score = float(cells.size()) / prev_cell_count;

    val = alpha * cell_increase_score + (1 - alpha) * size_score;
  }

  // do stuff here to update for the next scoring cycle
  prev_cell_count = size();

  return val;
}

const std::vector<int>& Bin::get_cells() const
{
  return cells;
}

void Bin::copy_untested_cells(
  const std::vector<int>& possible_cells, std::vector<int>& untested_cells)
{
  untested_cells.clear();

  lock_bin();
  for (int cell : possible_cells) {
    if (!std::binary_search(
          cells.begin(), cells.begin() + num_unique_cells, cell)) {
      untested_cells.push_back(cell);
    }
  }
  unlock_bin();
}
void Bin::lock_bin()
{
  omp_set_lock(&lock);
}

void Bin::unlock_bin()
{
  omp_unset_lock(&lock);
}
}; // namespace openmc