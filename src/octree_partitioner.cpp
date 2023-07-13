#include "openmc/octree_partitioner.h"
#include "openmc/error.h"
#include "openmc/random_dist.h"
#include "openmc/timer.h"

#include <omp.h>

#include <algorithm>
#include <fstream>
#include <queue>
#include <stack>

#include <stdlib.h>

namespace openmc {

constexpr int32_t OMCP_CURRENT_VERSION[3] = {1, 0, 0};

constexpr double REFINEMENT_SEARCH_DENSITY = 0.125;
constexpr double REFINEMENT_TIMEOUT = 10.0;
constexpr int32_t REFINEMENT_GRID_RES = 128; // ideally should be a power of 2

constexpr int32_t INFORMATION_REFILLING_START_DEPTH_OFFSET = 1;

enum OMCPStructureType {
  Octree = 1,
  KdTree = 2,            // not yet implemented in OMCP
  ZPlanePartitioner = 3, // not yet implemented in OMCP
};

struct OctreeConstructionTask {
  OctreeUncompressedNode* node_;
  std::vector<CellPoint> points_;

  OctreeConstructionTask() = default;
  OctreeConstructionTask(
    OctreeUncompressedNode* n, const std::vector<CellPoint>& p)
    : node_(n), points_(p)
  {}
};

OctreeNode::OctreeNode() : data_(0) {}

constexpr uint32_t OCTREE_LEAF_FLAG = (1 << 31);
void OctreeNode::mark_as_leaf()
{
  data_ |= OCTREE_LEAF_FLAG;
}

bool OctreeNode::is_leaf() const
{
  return (data_ & OCTREE_LEAF_FLAG);
}

void OctreeNode::store_data(uint32_t data)
{
  // remove everything except the flag
  this->data_ = (this->data_ & OCTREE_LEAF_FLAG) | data;
}

uint32_t OctreeNode::read_data() const
{
  return (data_ & ~OCTREE_LEAF_FLAG);
}

OctreeUncompressedNode::OctreeUncompressedNode()
  : children_(nullptr), parent_(nullptr), depth_(0)
{}

bool OctreeUncompressedNode::is_leaf() const
{
  return (children_ == nullptr);
}

void OctreeUncompressedNode::subdivide()
{
  AABB resultant_boxes[8];
  resultant_boxes[0] = box_;
  for (int i = 0; i < 3; i++) {
    AABB temp_box_buffer[8];

    int next_index = 0;
    for (int idx = 0; idx < (1 << i); idx++) {
      // split on i-th axis
      double midpoint = box_.get_center()[i];

      int j = ((i + 1) % 3);
      int k = ((i + 2) % 3);

      AABB splitted_boxes[2] {resultant_boxes[idx], resultant_boxes[idx]};

      splitted_boxes[0].max_[i] = midpoint;
      splitted_boxes[1].min_[i] = midpoint;

      temp_box_buffer[next_index++] = splitted_boxes[1];
      temp_box_buffer[next_index++] = splitted_boxes[0];
    }

    // move the results to the next splitting stage
    std::copy(temp_box_buffer, temp_box_buffer + (2 << i), resultant_boxes);
  }

  for (int i = 0; i < 8; i++) {
    children_[i].box_ = resultant_boxes[i];
    children_[i].depth_ = depth_ + 1;
    children_[i].parent_ = this;
  }
}

bool OctreeUncompressedNode::contains(int cell) const
{
  return std::binary_search(cells_.begin(), cells_.end(), cell);
}

void pick_untested_cells(const std::vector<int>& tested_cells,
  const std::vector<int>& possible_cells, std::vector<int>& untested_cells)
{
  untested_cells.clear();
  int next_idx = 0;
  for (int cell : possible_cells) {
    if (next_idx < tested_cells.size() && tested_cells[next_idx] == cell) {
      next_idx++;
    } else {
      untested_cells.push_back(cell);
    }
  }
}

void refine_octree_random(const Universe& univ,
  const UniversePartitioner& fallback, const AABB& bounds,
  const std::vector<OctreeUncompressedNode*>& leaves)
{
  const int32_t num_threads = omp_get_max_threads();

  // generate the seeds
  std::vector<uint64_t[2]> rng_node_selec(num_threads);
  std::vector<uint64_t[3]> rng_pos(num_threads);
  for (int i = 0; i < num_threads; i++) {
    rng_node_selec[i][0] = i;
    rng_node_selec[i][1] = i + num_threads;
    for (int j = 0; j < 3; j++) {
      rng_pos[i][j] = i + num_threads * (j + 2);
    }
  }

  using ProbBinT = ProbabilityBin<OctreeUncompressedNode>;
  std::vector<ProbBinT> prob_bin_grid(
    REFINEMENT_GRID_RES * REFINEMENT_GRID_RES * REFINEMENT_GRID_RES);

  Position prob_bin_dim;
  for (int i = 0; i < 3; i++) {
    prob_bin_dim[i] = (bounds.max_[i] - bounds.min_[i]) / REFINEMENT_GRID_RES;
  }

  for (auto leaf : leaves) {
    int idx = 0;
    for (int i = 0; i < 3; i++) {
      idx =
        REFINEMENT_GRID_RES * idx +
        int((leaf->box_.get_center()[i] - bounds.min_[i]) / prob_bin_dim[i]);
    }

    omp_lock_t dummy;
    prob_bin_grid[idx].contained_nodes_.emplace_back(leaf, dummy);
    for (int cell : leaf->cells_) {
      prob_bin_grid[idx].found_cells_.push_back(cell);
    }
  }

  for (auto& prob_bin : prob_bin_grid) {
    prob_bin.update_common_cells();
    for (auto& p : prob_bin.contained_nodes_) {
      omp_init_lock(&p.second);
    }
  }

  std::vector<double> bin_cdf(prob_bin_grid.size());

  Timer timeout_timer;
  timeout_timer.start();
  int iteration = 0;
  while (timeout_timer.elapsed() < REFINEMENT_TIMEOUT) {
    int num_search_points =
      static_cast<int>(bounds.volume() * REFINEMENT_SEARCH_DENSITY);
    int num_points_searched = 0;

    // first, generate cdf
    double total_cdf = 0.0;
    for (int i = 0; i < prob_bin_grid.size(); i++) {
      total_cdf += prob_bin_grid[i].compute_score();
      bin_cdf[i] = total_cdf;
    }

    for (double& cdf : bin_cdf) {
      cdf /= total_cdf;
    }

#pragma omp parallel for
    for (int tid = 0; tid < num_threads; tid++) {
      std::vector<int> untested_cells, temp_buf;
      untested_cells.reserve(256);
      temp_buf.reserve(256);

      // we don't need a lock for loop iteration logic since nothing bad will
      // happen if race conditions occur maybe we do a few points less or more,
      // but that doesn't matter
      while (num_points_searched < num_search_points) {
        num_points_searched++;

        ProbBinT* prob_bin;
        do {
          double cdf_val =
            uniform_distribution(0.0, 1.0, &rng_node_selec[tid][0]);
          auto cdf_iter =
            std::upper_bound(bin_cdf.begin(), bin_cdf.end(), cdf_val);
          size_t bin_idx = std::distance(bin_cdf.begin(), cdf_iter);
          if (bin_idx == bin_cdf.size()) {
            bin_idx--;
          }

          prob_bin = &prob_bin_grid[bin_idx];
        } while (prob_bin->contained_nodes_.size() == 0);
        prob_bin->num_searched_++;

        size_t idx = prob_bin->pick_node(&rng_node_selec[tid][1]);
        if (idx >= prob_bin->contained_nodes_.size()) {
          idx = prob_bin->contained_nodes_.size() - 1;
        }

        OctreeUncompressedNode* current = prob_bin->contained_nodes_[idx].first;
        omp_lock_t* cell_lock = &prob_bin->contained_nodes_[idx].second;

        CellPointUncompressed point;
        for (int i = 0; i < 3; i++) {
          point.pos_[i] = uniform_distribution(
            current->box_.min_[i], current->box_.max_[i], &rng_pos[tid][i]);
        }

        omp_set_lock(cell_lock);
        point.cell_ = univ.find_cell_for_point(current->cells_, point.pos_);
        omp_unset_lock(cell_lock);

        if (point.cell_ == -1) {
          prob_bin->num_found_++;

          omp_set_lock(cell_lock);
          pick_untested_cells(
            current->cells_, prob_bin->common_cells_, untested_cells);
          omp_unset_lock(cell_lock);

          point.cell_ = univ.find_cell_for_point(untested_cells, point.pos_);
          if (point.cell_ == -1) {
            point.cell_ =
              univ.find_cell_for_point(current->parent_->cells_, point.pos_);

            if (point.cell_ == -1) {
              Direction dummy {0, 0, 1};
              const auto& possible_cells =
                fallback.get_cells(point.pos_, dummy);

              pick_untested_cells(untested_cells, possible_cells, temp_buf);

              omp_set_lock(cell_lock);
              pick_untested_cells(current->cells_, temp_buf, untested_cells);
              omp_unset_lock(cell_lock);

              point.cell_ =
                univ.find_cell_for_point(untested_cells, point.pos_);

              // very rarely, even the fallback misses a cell
              // we need to do an exhaustive search or the program will segfault
              if (point.cell_ == -1) {
                point.cell_ = univ.find_cell_for_point(univ.cells_, point.pos_);
              }
            }

            prob_bin->add_cell(point.cell_);
          }

          omp_set_lock(cell_lock);
          // insertion sort
          current->cells_.push_back(point.cell_);
          for (int i = current->cells_.size() - 1; i > 0; i--) {
            if (current->cells_[i] < current->cells_[i - 1]) {
              std::swap(current->cells_[i], current->cells_[i - 1]);
            }
          }
          omp_unset_lock(cell_lock);
        }
      }
    }
    iteration++;
  }
}

OctreePartitioner::OctreePartitioner(
  const Universe& univ, int target_cells_per_node)
  : fallback(univ)
{
  constexpr double half_side_length = 130.0;
  bounds_.min_ =
    Position(-half_side_length, -half_side_length, -half_side_length);
  bounds_.max_ = Position(half_side_length, half_side_length, half_side_length);

  if (univ.cells_.size() <= target_cells_per_node) {
    warning("Universe has only " + std::to_string(univ.cells_.size()) +
            " cells, which is below the target cells per node, which is " +
            std::to_string(target_cells_per_node) +
            " cells. Octree will only consist of root.");

    OctreeNode node;
    node.store_data(0);
    node.mark_as_leaf();
    cell_data_.push_back(univ.cells_);
    return;
  }

  write_message("Building octree...", 5);

  Timer construction_timer;
  construction_timer.start();

  auto points_in_bounds =
    binned_point_search<CellPoint>(univ, fallback, bounds_);

  OctreeUncompressedNode root;
  root.box_ = bounds_;
  root.depth_ = 0;

  std::sort(points_in_bounds.begin(), points_in_bounds.end());
  int prev_cell = -1;
  for (const auto& p : points_in_bounds) {
    if (prev_cell != p.get_cell()) {
      root.cells_.push_back(p.get_cell());
      prev_cell = p.get_cell();
    }
  }
  root.num_unique_cells_ = root.cells_.size();

  num_nodes_ = 1;
  num_leaves_ = 0;

  double depth_vh_mult[] = {1.0, 1.0, 1.0, 1.5, 2.5, 4.0, 6.0, 12.5, 19.0, 32.0,
    64.0, 128.0, 999.0, 9999.0, 99999.0};

  NodeAllocator<OctreeUncompressedNode, 8> node_alloc;
  std::vector<OctreeUncompressedNode*> nodes_to_propagate {&root};
  std::vector<OctreeUncompressedNode*> leaves;

  // this section of code still needs to be multithreaded
  // it can become a bottleneck, espcially with large number of points
  std::queue<OctreeConstructionTask> unprocessed_tasks;
  unprocessed_tasks.emplace(&root, points_in_bounds);
  while (!unprocessed_tasks.empty()) {

    auto cur_task = std::move(unprocessed_tasks.front());
    unprocessed_tasks.pop();

    // subdivide
    cur_task.node_->children_ = node_alloc.allocate();
    cur_task.node_->subdivide();

    OctreeConstructionTask child_tasks[8];
    for (int i = 0; i < 8; i++) {
      child_tasks[i].node_ = &cur_task.node_->children_[i];
      nodes_to_propagate.push_back(child_tasks[i].node_);
    }

    // sort points
    for (const auto& point : cur_task.points_) {
      child_tasks[point.get_child_index(cur_task.node_->depth_)]
        .points_.push_back(point);
    }

    double parent_vh =
      cur_task.node_->box_.volume() * cur_task.node_->num_unique_cells_;
    double children_vh = 0.0;

    // post processing (make nodes leaves or push on construction stack)
    bool force_subdiv = false;
    for (int i = 0; i < 8; i++) {
      // count the number of unique cells
      int num_unique_cells = 0;
      int prev_cell = -1;
      for (const auto& p : child_tasks[i].points_) {
        if (p.get_cell() != prev_cell) {
          num_unique_cells++;
          prev_cell = p.get_cell();
        }
      }

      child_tasks[i].node_->num_unique_cells_ = num_unique_cells;

      children_vh += child_tasks[i].node_->box_.volume() * num_unique_cells;
      if (num_unique_cells > target_cells_per_node) {
        force_subdiv = true;
      }
    }

    if (force_subdiv ||
        depth_vh_mult[cur_task.node_->depth_] * children_vh < parent_vh) {
      // continue subdivision on this branch
      num_nodes_ += 8;
      for (int i = 0; i < 8; i++) {
        child_tasks[i].node_->cells_.reserve(
          child_tasks[i].node_->num_unique_cells_);
        prev_cell = -1;
        for (const auto& p : child_tasks[i].points_) {
          if (p.get_cell() != prev_cell) {
            child_tasks[i].node_->cells_.push_back(p.get_cell());
            prev_cell = p.get_cell();
          }
        }

        unprocessed_tasks.push(std::move(child_tasks[i]));
      }
    } else {
      // terminate subdivision on this branch
      cur_task.node_->children_ = nullptr;
      leaves.push_back(cur_task.node_);
      num_leaves_++;
    }

    // free memory
    cur_task.points_.clear();
    cur_task.points_.shrink_to_fit();
  }

  refine_octree_random(univ, fallback, bounds_, leaves);

  // now, build the cells list
  struct NodeComp {
    inline bool operator()(const OctreeUncompressedNode* lhs,
      const OctreeUncompressedNode* rhs) const
    {
      return (lhs->depth_ < rhs->depth_);
    }
  };
  std::sort(nodes_to_propagate.rbegin(), nodes_to_propagate.rend(), NodeComp());

  for (auto ptr : nodes_to_propagate) {
    auto& cur = *ptr;

    if (!cur.is_leaf()) {
      std::sort(cur.cells_.begin(), cur.cells_.end());

      int prev_cell = cur.cells_[0];

      int next_idx = 1;
      for (int i = 1; i < cur.cells_.size(); i++) {
        if (cur.cells_[i] != prev_cell) {
          cur.cells_[next_idx] = cur.cells_[i];
          next_idx++;
          prev_cell = cur.cells_[i];
        }
      }
      cur.cells_.resize(next_idx);
      cur.cells_.shrink_to_fit();
    }

    if (cur.parent_) {
      for (int cell : cur.cells_) {
        cur.parent_->cells_.push_back(cell);
      }
    }
  }

  // a possible way to remove issues with missing information in nodes
  std::stack<OctreeUncompressedNode*> uninfo_refilled_nodes;
  uninfo_refilled_nodes.push(&root);
  while (!uninfo_refilled_nodes.empty()) {
    auto cur = uninfo_refilled_nodes.top();
    uninfo_refilled_nodes.pop();

    bool should_collect_leaves = false;
    for (int i = 0; i < 8; i++) {
      if (cur->children_[i].is_leaf()) {
        should_collect_leaves = true;
        break;
      }
    }

    auto collection_start = cur;
    for (int i = 0; i < INFORMATION_REFILLING_START_DEPTH_OFFSET; i++) {
      if (collection_start->parent_) {
        collection_start = collection_start->parent_;
      }
    }

    const auto& unique_cells = collection_start->cells_;

    for (int i = 0; i < 8; i++) {
      if (cur->children_[i].is_leaf()) {
        int orig_size = cur->children_[i].cells_.size();

        auto search_start = cur->children_[i].cells_.begin();
        auto search_end = cur->children_[i].cells_.begin() + orig_size;

        cur->children_[i].num_original_cells_ = orig_size;

        for (int cell : unique_cells) {
          if (!std::binary_search(search_start, search_end, cell)) {
            cur->children_[i].cells_.push_back(cell);
          }
        }
      } else {
        uninfo_refilled_nodes.push(&cur->children_[i]);
      }
    }
  }

  // now copy everything to array
  nodes_.reserve(num_nodes_);
  cell_data_.reserve(num_leaves_);
  orig_size_.reserve(num_leaves_);
  std::queue<const OctreeUncompressedNode*> unwritten_nodes;
  unwritten_nodes.push(&root);
  while (!unwritten_nodes.empty()) {
    auto cur = unwritten_nodes.front();
    unwritten_nodes.pop();

    // convert and write
    OctreeNode compressed;

    if (cur->is_leaf()) {
      compressed.store_data(cell_data_.size());
      cell_data_.push_back(std::move(cur->cells_));
      compressed.mark_as_leaf();

      orig_size_.push_back(cur->num_original_cells_);
    } else {
      compressed.store_data(nodes_.size() + 1 + unwritten_nodes.size());
      for (int i = 0; i < 8; i++) {
        unwritten_nodes.push(&cur->children_[i]);
      }
    }

    nodes_.push_back(compressed);
  }

  write_message("Octree construction completed in " +
                  std::to_string(construction_timer.elapsed()) + " seconds.",
    5);
}

OctreePartitioner::OctreePartitioner(
  const Universe& univ, const std::string& path)
  : fallback(univ)
{
  write_message("Reading octree from " + path, 5);

  std::fstream octree_file(path, std::ios::in | std::ios::binary);
#define READ_BINARY(x) octree_file.read((char*)&x, sizeof(x))

  for (int i = 0; i < 3; i++) {
    int version;
    READ_BINARY(version);
    if (version != OMCP_CURRENT_VERSION[i]) {
      fatal_error(
        "OpenMC cannot read an unsupported OMCP file version! Please note that "
        "OpenMC currently cannot read older file versions!");
    }
  }

  OMCPStructureType structure_type;
  READ_BINARY(structure_type);

  if (structure_type != OMCPStructureType::Octree) {
    fatal_error(
      "OpenMC currently only supports octrees for the OMCP file format!");
  }

  std::vector<uint16_t> comp_cell_data;

  READ_BINARY(bounds_);

  READ_BINARY(num_nodes_);

  int num_leaves = 0;
  nodes_.resize(num_nodes_);
  for (auto& node : nodes_) {
    READ_BINARY(node);
    num_leaves++;
  }
  cell_data_.reserve(num_leaves);

  int num_cells;
  READ_BINARY(num_cells);
  comp_cell_data.resize(num_cells);
  for (auto& cell : comp_cell_data) {
    READ_BINARY(cell);
  }

  int next_cell_index = 0;
  for (auto& node : nodes_) {
    if (node.is_leaf()) {
      std::vector<int> cells;
      cells.resize(node.read_data());

      for (int& cell : cells) {
        cell = comp_cell_data.at(next_cell_index);
        next_cell_index++;
      }

      node.store_data(cell_data_.size());
      cell_data_.push_back(std::move(cells));
    }
  }

  refill_information();
}

OctreePartitioner::~OctreePartitioner() {}

const vector<int32_t>& OctreePartitioner::get_cells(
  Position r, Direction u) const
{
  if (!bounds_.contains(r)) {
    // discount this get_cells call from the stats to only focus on points
    // within the octree
    return get_cells_fallback(r, u);
  }

  Position node_dim;
  for (int i = 0; i < 3; i++) {
    node_dim[i] = (bounds_.max_[i] - bounds_.min_[i]) * 0.5;
  }

  Position center = bounds_.get_center();
  auto current = nodes_[0];
  while (!current.is_leaf()) {
    // halve the node dim
    int idx = 0;
    for (int i = 0; i < 3; i++) {
      node_dim[i] *= 0.5;
      bool less = (r[i] < center[i]);
      center[i] += node_dim[i] * (less ? -1 : 1);
      idx = 2 * idx + int(less);
    }

    current = nodes_[current.read_data() + idx];
  }

  const auto& cells = cell_data_[current.read_data()];
  if (cells.empty()) {
    return get_cells_fallback(r, u);
  }

  return cells;
}

const vector<int32_t>& OctreePartitioner::get_cells_fallback(
  Position r, Direction u) const
{
  return fallback.get_cells(r, u);
}

void OctreePartitioner::export_to_file(const std::string& path) const
{
  std::vector<uint16_t> comp_cell_data;
  for (int i = 0; i < cell_data_.size(); i++) {
    for (int j = 0; j < orig_size_[i]; j++) {
      comp_cell_data.push_back(cell_data_[i][j]);
    }
  }

  std::fstream octree_file(path, std::ios::out | std::ios::binary);

#define WRITE_BINARY(x) octree_file.write((char*)&x, sizeof(x))

  for (int i = 0; i < 3; i++) {
    WRITE_BINARY(OMCP_CURRENT_VERSION[i]);
  }
  OMCPStructureType structure_type = OMCPStructureType::Octree;
  WRITE_BINARY(structure_type);

  WRITE_BINARY(bounds_);

  WRITE_BINARY(num_nodes_);
  for (const auto& raw_node : nodes_) {
    auto node = raw_node;
    if (node.is_leaf()) {
      node.store_data(orig_size_[node.read_data()]);
    }
    WRITE_BINARY(node);
  }

  int num_cells = comp_cell_data.size();
  WRITE_BINARY(num_cells);
  for (auto cell : comp_cell_data) {
    WRITE_BINARY(cell);
  }

  write_message("Exported octree to " + path, 5);
}

// this method works on the compressed octree, not hte uncompressed one
// since it is pretty slow (and for some reason creates an octree that has the
// same failure rate but is slower), only use it if you are reading from a file
void OctreePartitioner::refill_information()
{
  struct OctreeNodeExtraInfo {
    OctreeNodeExtraInfo()
      : parent_(nullptr), eiparent_(nullptr), current_(nullptr),
        eichildren_(nullptr), children_(nullptr), depth_(0)
    {}

    OctreeNode* parent_;
    OctreeNode* current_;
    OctreeNode* children_;

    OctreeNodeExtraInfo* eiparent_;
    OctreeNodeExtraInfo* eichildren_;

    std::vector<int> cells_;

    uint32_t depth_ = 0;

    bool is_leaf() const { return (children_ != nullptr); }
  };

  // make our nodes easier to work with
  std::vector<OctreeNodeExtraInfo> einodes(nodes_.size());
  std::vector<OctreeNodeExtraInfo*> nodes_to_propagate(nodes_.size());
  for (int i = 0; i < nodes_.size(); i++) {
    OctreeNodeExtraInfo* einode = &einodes[i];
    nodes_to_propagate[i] = einode;

    einode->current_ = &nodes_[i];

    if (!nodes_[i].is_leaf()) {
      int idx = nodes_[i].read_data();
      einode->children_ = &nodes_[idx];
      einode->eichildren_ = &einodes[idx];

      for (int i = 0; i < 8; i++) {
        einode->eichildren_[i].parent_ = einode->current_;
        einode->eichildren_[i].eiparent_ = einode;
        einode->eichildren_[i].depth_ = einode->depth_ + 1;
      }
    } else {
      einode->cells_ = std::move(cell_data_[nodes_[i].read_data()]);
    }
  }

  // propagate all cells forward
  struct DepthComp {
    inline bool operator()(
      const OctreeNodeExtraInfo* lhs, const OctreeNodeExtraInfo* rhs) const
    {
      return (lhs->depth_ < rhs->depth_);
    }
  };
  std::sort(
    nodes_to_propagate.rbegin(), nodes_to_propagate.rend(), DepthComp());
  for (auto ptr : nodes_to_propagate) {
    auto& cur = *ptr;
    if (!cur.is_leaf()) {
      std::sort(cur.cells_.begin(), cur.cells_.end());

      int prev_cell = cur.cells_[0];

      int next_idx = 1;
      for (int i = 1; i < cur.cells_.size(); i++) {
        if (cur.cells_[i] != prev_cell) {
          cur.cells_[next_idx] = cur.cells_[i];
          next_idx++;
          prev_cell = cur.cells_[i];
        }
      }
      cur.cells_.resize(next_idx);
      cur.cells_.shrink_to_fit();
    }

    if (cur.parent_) {
      for (int cell : cur.cells_) {
        cur.eiparent_->cells_.push_back(cell);
      }
    }
  }

  // now propagate all cells downward
  std::stack<OctreeNodeExtraInfo*> uninfo_refilled_nodes;
  uninfo_refilled_nodes.push(nodes_to_propagate.back());
  while (!uninfo_refilled_nodes.empty()) {
    auto cur = uninfo_refilled_nodes.top();
    uninfo_refilled_nodes.pop();

    bool should_collect_leaves = false;
    for (int i = 0; i < 8; i++) {
      if (cur->children_[i].is_leaf()) {
        should_collect_leaves = true;
        break;
      }
    }

    auto collection_start = cur;
    for (int i = 0; i < INFORMATION_REFILLING_START_DEPTH_OFFSET; i++) {
      if (collection_start->parent_) {
        collection_start = collection_start->eiparent_;
      }
    }

    const auto& unique_cells = collection_start->cells_;

    for (int i = 0; i < 8; i++) {
      if (cur->children_[i].is_leaf()) {
        int orig_size = cur->eichildren_[i].cells_.size();

        auto search_start = cur->eichildren_[i].cells_.begin();
        auto search_end = cur->eichildren_[i].cells_.begin() + orig_size;

        for (int cell : unique_cells) {
          if (!std::binary_search(search_start, search_end, cell)) {
            cur->eichildren_[i].cells_.push_back(cell);
          }
        }
        // now that downpropagation is done for this node, update cell_data
        cell_data_[cur->eichildren_[i].current_->read_data()] =
          std::move(cur->eichildren_[i].cells_);
      } else {
        uninfo_refilled_nodes.push(&cur->eichildren_[i]);
      }
    }
  }
}

} // namespace openmc