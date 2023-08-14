#include "openmc/octree_partitioner.h"
#include "openmc/error.h"
#include "openmc/partitioners.h"
#include "openmc/random_dist.h"
#include "openmc/timer.h"

#include <omp.h>

#include <algorithm>
#include <fstream>
#include <queue>
#include <stack>

#include <stdlib.h>

namespace openmc {
constexpr double REFINEMENT_SEARCH_DENSITY = 0.125;
constexpr double REFINEMENT_TIMEOUT = 10.0;
constexpr int32_t REFINEMENT_GRID_RES = 128; // ideally should be a power of 2

constexpr int32_t INFORMATION_REFILLING_START_DEPTH_OFFSET = 1;

// This structure is used during construction to assocaite nodes with arrays of
// points. I dont directly include the points_ within the nodes because it
// appears to lower performance. Perhaps that is due to the reduced cache hit
// rate during refinement.
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

// 32nd bit is on, all else is off
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
  // Store the first 31 bits while keeping the value of the leaf flag
  this->data_ = (this->data_ & OCTREE_LEAF_FLAG) | data;
}

uint32_t OctreeNode::read_data() const
{
  // Return all bits except the leaf flag
  return (data_ & ~OCTREE_LEAF_FLAG);
}

void OctreeNode::store_raw_data(uint32_t val)
{
  data_ = val;
}

uint32_t OctreeNode::read_raw_data()
{
  return data_;
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
  // Recursively splits list of boxes on each axis to get 8 sub boxes
  // Initially this list only consists of parent's box

  BoundingBox resultant_boxes[8];
  resultant_boxes[0] = box_;
  for (int i = 0; i < 3; i++) {
    BoundingBox temp_box_buffer[8];

    // Split each box in resultant_boxes on axis i and store them in
    // temp_box_buffer
    int next_index = 0;
    for (int idx = 0; idx < (1 << i); idx++) {
      // split on i-th axis
      double midpoint = box_.get_center()[i];

      int j = ((i + 1) % 3);
      int k = ((i + 2) % 3);

      BoundingBox splitted_boxes[2] {resultant_boxes[idx], resultant_boxes[idx]};

      splitted_boxes[0].max_[i] = midpoint;
      splitted_boxes[1].min_[i] = midpoint;

      temp_box_buffer[next_index++] = splitted_boxes[1];
      temp_box_buffer[next_index++] = splitted_boxes[0];
    }

    // move the results to the next splitting stage
    std::copy(temp_box_buffer, temp_box_buffer + (2 << i), resultant_boxes);
  }

  // Move boxes over to children
  for (int i = 0; i < 8; i++) {
    children_[i].box_ = resultant_boxes[i];
    children_[i].depth_ = depth_ + 1;
    children_[i].parent_ = this;
  }
}

// Assumes that cells are storted
bool OctreeUncompressedNode::contains(int cell) const
{
  return std::binary_search(cells_.begin(), cells_.end(), cell);
}

// Serves same functionality as Bin::copy_untested_cells. See the comment on
// that in partitioner_utils.h
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

// The refinement process picks a random node within the octree and samples a
// point instead it. If the leaf does not have the cell, then it utilizes the
// fallback to find the actual cell and then inserts the new cell into the leaf.
// Note that it does not subdivide a leaf if the number of cells exceeds the
// target cells per leaf paramter. This is done in order to prevent pockets of
// missing information from forming.
void refine_octree_random(const Universe& univ,
  const UniversePartitioner& fallback, const BoundingBox& bounds,
  const std::vector<OctreeUncompressedNode*>& leaves)
{
  // Generate the seeds that each thread will use ahead of time
  const int32_t num_threads = omp_get_max_threads();
  std::vector<uint64_t[2]> rng_node_selec(num_threads);
  std::vector<uint64_t[3]> rng_pos(num_threads);
  for (int i = 0; i < num_threads; i++) {
    rng_node_selec[i][0] = i;
    rng_node_selec[i][1] = i + num_threads;
    for (int j = 0; j < 3; j++) {
      rng_pos[i][j] = i + num_threads * (j + 2);
    }
  }

  // Create our probability bins. These basically help us importance sample
  // nodes for refinement.
  using ProbBinT = ProbabilityBin<OctreeUncompressedNode>;
  std::vector<ProbBinT> prob_bin_grid(
    REFINEMENT_GRID_RES * REFINEMENT_GRID_RES * REFINEMENT_GRID_RES);
  Position prob_bin_dim;
  for (int i = 0; i < 3; i++) {
    prob_bin_dim[i] = (bounds.max_[i] - bounds.min_[i]) / REFINEMENT_GRID_RES;
  }

  // For each leaf, add it to a probabilty bin
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

  // For each bin, initalize the locks
  for (auto& prob_bin : prob_bin_grid) {
    prob_bin.update_common_cells();
    for (auto& p : prob_bin.contained_nodes_) {
      omp_init_lock(&p.second);
    }
  }

  Timer timeout_timer;
  timeout_timer.start();

  std::vector<double> bin_cmf(prob_bin_grid.size());
  int iteration = 0;

  // We terminate the loop after refining for a certain amount of time
  while (timeout_timer.elapsed() < REFINEMENT_TIMEOUT) {
    int num_search_points =
      static_cast<int>(bounds.volume() * REFINEMENT_SEARCH_DENSITY);
    int num_points_searched = 0;

    // First, generate cmf
    double total_cmf = 0.0;
    for (int i = 0; i < prob_bin_grid.size(); i++) {
      total_cmf += prob_bin_grid[i].compute_score();
      bin_cmf[i] = total_cmf;
    }

    for (double& cmf : bin_cmf) {
      cmf /= total_cmf;
    }

    // Launch a job for each thread
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

        // Find a probabity bin
        ProbBinT* prob_bin;
        do {
          double cmf_val =
            uniform_distribution(0.0, 1.0, &rng_node_selec[tid][0]);
          auto cmf_iter =
            std::upper_bound(bin_cmf.begin(), bin_cmf.end(), cmf_val);
          size_t bin_idx = std::distance(bin_cmf.begin(), cmf_iter);
          if (bin_idx == bin_cmf.size()) {
            bin_idx--;
          }

          prob_bin = &prob_bin_grid[bin_idx];
        } while (prob_bin->contained_nodes_.size() == 0);
        prob_bin->num_searched_++;

        // Now pick a node
        size_t idx = prob_bin->pick_node(&rng_node_selec[tid][1]);
        if (idx >= prob_bin->contained_nodes_.size()) {
          idx = prob_bin->contained_nodes_.size() - 1;
        }

        OctreeUncompressedNode* current = prob_bin->contained_nodes_[idx].first;
        omp_lock_t* cell_lock = &prob_bin->contained_nodes_[idx].second;

        // Generate a point within in
        CellPointUncompressed point;
        for (int i = 0; i < 3; i++) {
          point.pos_[i] = uniform_distribution(
            current->box_.min_[i], current->box_.max_[i], &rng_pos[tid][i]);
        }

        // Search the point for the node
        omp_set_lock(cell_lock);
        point.cell_ = univ.find_cell_for_point(current->cells_, point.pos_);
        omp_unset_lock(cell_lock);

        if (point.cell_ == -1) {
          // Point was not found, let's search in the bin's list of common cells
          prob_bin->num_found_++;

          omp_set_lock(cell_lock);
          pick_untested_cells(
            current->cells_, prob_bin->common_cells_, untested_cells);
          omp_unset_lock(cell_lock);

          point.cell_ = univ.find_cell_for_point(untested_cells, point.pos_);
          if (point.cell_ == -1) {
            // Still not found, let's check the parent
            point.cell_ =
              univ.find_cell_for_point(current->parent_->cells_, point.pos_);

            if (point.cell_ == -1) {
              // We are going to have to check the fallback's list
              Direction dummy {0, 0, 1};
              const auto& possible_cells =
                fallback.get_cells(point.pos_, dummy);

              pick_untested_cells(untested_cells, possible_cells, temp_buf);

              omp_set_lock(cell_lock);
              pick_untested_cells(current->cells_, temp_buf, untested_cells);
              omp_unset_lock(cell_lock);

              point.cell_ =
                univ.find_cell_for_point(untested_cells, point.pos_);

              if (point.cell_ == -1) {
                // very rarely, even the fallback misses a cell
                // we need to do an exhaustive search or the program will
                // segfault
                point.cell_ = univ.find_cell_for_point(univ.cells_, point.pos_);
              }
            }

            // Add the missing cell to the bin
            prob_bin->add_cell(point.cell_);
          }

          // insertion sort into the leaf's list of cells
          omp_set_lock(cell_lock);
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
  const Universe& univ, const BoundingBox& bounds, int target_cells_per_node)
  : fallback(univ), bounds_(bounds)
{

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

  // Search for points
  auto points_in_bounds =
    binned_point_search<CellPoint>(univ, fallback, bounds_);

  // Set up our root
  OctreeUncompressedNode root;
  root.box_ = bounds_;
  root.depth_ = 0;

  // Set up the number of unique cells here.
  // We need to sort anyway for construction.
  std::sort(points_in_bounds.begin(), points_in_bounds.end());
  int prev_cell = -1;
  for (const auto& p : points_in_bounds) {
    if (prev_cell != p.get_cell()) {
      root.cells_.push_back(p.get_cell());
      prev_cell = p.get_cell();
    }
  }
  root.num_unique_cells_ = root.cells_.size();

  // The root already is the first node
  num_nodes_ = 1;
  num_leaves_ = 0;

  // Paramters to control subdivision
  double depth_vh_mult[] = {1.0, 1.0, 1.0, 1.5, 2.5, 4.0, 6.0, 12.5, 19.0, 32.0,
    64.0, 128.0, 999.0, 9999.0, 99999.0};

  // Node allocator to manage uncompressed child node allocation
  NodeAllocator<OctreeUncompressedNode, 8> node_alloc;
  // A list of all leaves, used later during refinement
  std::vector<OctreeUncompressedNode*> leaves;
  // List of all nodes, used later during information refilling
  std::vector<OctreeUncompressedNode*> nodes_to_propagate {&root};

  // this section of code still needs to be multithreaded.
  // it can become a bottleneck, especially with large number of points.
  std::queue<OctreeConstructionTask> unprocessed_tasks;
  unprocessed_tasks.emplace(&root, points_in_bounds);
  while (!unprocessed_tasks.empty()) {

    auto cur_task = std::move(unprocessed_tasks.front());
    unprocessed_tasks.pop();

    // allocate and subdivide subdivide
    cur_task.node_->children_ = node_alloc.allocate();
    cur_task.node_->subdivide();

    OctreeConstructionTask child_tasks[8];
    for (int i = 0; i < 8; i++) {
      child_tasks[i].node_ = &cur_task.node_->children_[i];
      nodes_to_propagate.push_back(child_tasks[i].node_);
    }

    // Sort points
    for (const auto& point : cur_task.points_) {
      child_tasks[point.get_child_index(cur_task.node_->depth_)]
        .points_.push_back(point);
    }

    // Initialize the VH values
    double parent_vh =
      cur_task.node_->box_.volume() * cur_task.node_->num_unique_cells_;
    double children_vh = 0.0;

    // determine the VH of a a split or determine if we have to split anyway
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
        // If we have too many cells, we have to split anyway
        force_subdiv = true;
      }
    }

    if (force_subdiv ||
        depth_vh_mult[cur_task.node_->depth_] * children_vh < parent_vh) {
      // Continue subdivision on this branch
      num_nodes_ += 8;
      for (int i = 0; i < 8; i++) {
        // Build the cell list and push to construction stack
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
      // Terminate subdivision on this branch
      cur_task.node_->children_ = nullptr;
      leaves.push_back(cur_task.node_);
      num_leaves_++;
    }

    // Free memory
    cur_task.points_.clear();
    cur_task.points_.shrink_to_fit();
  }

  // Run refinement
  refine_octree_random(univ, fallback, bounds_, leaves);

  // The next stage is information refilling. The first stage is forward
  // propagation, where we use the list of unique cells in the leaves to build
  // the list of unique cells for all nodes in the tree. The second stage is
  // down propagation, which takes those list of unique cells in the upper
  // levels of the tree (i.e. non-leaf nodes) and inserts them into the lower
  // levels of the tree.

  // First, clear out the cells in any upper levels of tree. These contain
  // information that is no longer useful since we now have refined information
  // in the lower levels of the tree. If we are not refining, we can entirely
  // skip the first stage of information refilling as the actual construction
  // already creates cell lists in the upper levels of the tree.
  for (auto ptr : nodes_to_propagate) {
    if (!ptr->is_leaf()) {
      ptr->cells_.clear();
    }
  }

  // Reverse sort the list of nodes according to their depth. This allows us to
  // process child nodes before thier parents nodes later.
  struct NodeComp {
    inline bool operator()(const OctreeUncompressedNode* lhs,
      const OctreeUncompressedNode* rhs) const
    {
      return (lhs->depth_ < rhs->depth_);
    }
  };
  std::sort(nodes_to_propagate.rbegin(), nodes_to_propagate.rend(), NodeComp());
  // Process nodes by depth.
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

  // Now down propagate the cell lists. We traverse through the tree and if we
  // come to an upper level node that has a leaf child, we down propagate the
  // cells. We take the upper level node's list of cells (or one of its
  // ancestors) and insert it into the cells of the tree.
  std::stack<OctreeUncompressedNode*> uninfo_refilled_nodes;
  uninfo_refilled_nodes.push(&root);
  while (!uninfo_refilled_nodes.empty()) {
    auto cur = uninfo_refilled_nodes.top();
    uninfo_refilled_nodes.pop();

    // Check if there is any leaf child
    bool should_down_propagate = false;
    for (int i = 0; i < 8; i++) {
      if (cur->children_[i].is_leaf()) {
        should_down_propagate = true;
        break;
      }
    }

    // Go up a few levels for the collection point. We might do this because the
    // parent itself might suffer from pockets of missing information.
    auto collection_point = cur;
    for (int i = 0; i < INFORMATION_REFILLING_START_DEPTH_OFFSET; i++) {
      if (collection_point->parent_) {
        collection_point = collection_point->parent_;
      }
    }

    for (int i = 0; i < 8; i++) {
      if (cur->children_[i].is_leaf()) {
        // Insert the list of nodes into this leaf

        cur->children_[i].num_original_cells_ = cur->children_[i].cells_.size();
        for (int cell : collection_point->cells_) {
          // Since the leaf's list of cells is sorted, we can binary search it
          // to avoid adding redundant cells
          if (!std::binary_search(cur->children_[i].cells_.begin(),
                cur->children_[i].cells_.begin() +
                  cur->children_[i].num_original_cells_,
                cell)) {
            cur->children_[i].cells_.push_back(cell);
          }
        }
      } else {
        // Push it onto the stack
        uninfo_refilled_nodes.push(&cur->children_[i]);
      }
    }
  }

  // Now reduce the information stored in each node as we no longer need it
  nodes_.reserve(num_nodes_);
  cell_data_.reserve(num_leaves_);
  // This is used to only output the non-information refilled cells to the file
  orig_size_.reserve(num_leaves_);

  // Go through each node and reduce the information
  std::queue<const OctreeUncompressedNode*> unwritten_nodes;
  unwritten_nodes.push(&root);
  while (!unwritten_nodes.empty()) {
    auto cur = unwritten_nodes.front();
    unwritten_nodes.pop();

    OctreeNode compressed;
    if (cur->is_leaf()) {
      // Store what index in cell_data_ the current leaf's cell list will be
      // located
      compressed.store_data(cell_data_.size());
      cell_data_.push_back(std::move(cur->cells_));
      compressed.mark_as_leaf();
      // Reduce size when exporting to file
      orig_size_.push_back(cur->num_original_cells_);
    } else {
      // Determine where the first child will be. We can only do this is we use
      // a queue instead of a stack.
      compressed.store_data(nodes_.size() + 1 + unwritten_nodes.size());
      for (int i = 0; i < 8; i++) {
        // Compress this node's children nodes too.
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
  const Universe& univ, const BoundingBox& bounds, hid_t file)
  : fallback(univ), bounds_(bounds)
{

  // read the nodes
  hid_t nodes_group = open_group(file, "node_data");
  read_attr_int(nodes_group, "n_nodes", &num_nodes_);
  std::vector<uint32_t> serialized_nodes;
  read_dataset(nodes_group, "nodes", serialized_nodes);
  close_group(nodes_group);

  // read the cells
  hid_t cells_group = open_group(file, "cell_data");
  int num_cells;
  read_attr_int(cells_group, "n_cells", &num_cells);
  std::vector<int32_t> comp_cell_data;
  read_dataset(cells_group, "cells", comp_cell_data);
  close_group(cells_group);

  // now take the information we read from the files and init our octree
  num_leaves_ = 0;
  nodes_.resize(num_nodes_);

  int next_cell_index = 0;

  for (int i = 0; i < num_nodes_; i++) {
    nodes_[i].store_raw_data(serialized_nodes[i]);
    if (nodes_[i].is_leaf()) {
      num_leaves_++;

      std::vector<int> cells;
      cells.resize(nodes_[i].read_data());

      for (int& cell : cells) {
        cell = comp_cell_data.at(next_cell_index);
        next_cell_index++;
      }

      nodes_[i].store_data(cell_data_.size());
      cell_data_.push_back(std::move(cells));
    }
  }

  refill_information();
}

OctreePartitioner::~OctreePartitioner() {}

void OctreePartitioner::export_to_hdf5(const std::string& path) const
{
  std::vector<int32_t> comp_cell_data;
  for (int i = 0; i < cell_data_.size(); i++) {
    for (int j = 0; j < orig_size_[i]; j++) {
      comp_cell_data.push_back(cell_data_[i][j]);
    }
  }

  hid_t file = file_open(path, 'w');

  // Write header
  write_attribute(file, "filetype", "partitioner");

  // write general partitioner information
  write_attribute(
    file, "part_type", static_cast<int>(PartitionerTypeID::Octree));
  write_dataset(file, "bounds_max", bounds_.max_);
  write_dataset(file, "bounds_min", bounds_.min_);

  // write the nodes
  hid_t nodes_group = create_group(file, "node_data");
  write_attribute(nodes_group, "n_nodes", nodes_.size());
  // Serialize the data
  std::vector<uint32_t> serialized_nodes(nodes_.size());
  for (int i = 0; i < nodes_.size(); i++) {
    auto node = nodes_[i]; // copy
    if (node.is_leaf()) {
      node.store_data(orig_size_[node.read_data()]);
    }
    serialized_nodes[i] = node.read_raw_data();
  }
  write_dataset(nodes_group, "nodes", serialized_nodes);
  close_group(nodes_group);

  // write the cell data
  hid_t cells_group = create_group(file, "cell_data");
  write_attribute(cells_group, "n_cells", comp_cell_data.size());
  write_dataset(cells_group, "cells", comp_cell_data);
  close_group(cells_group);

  file_close(file);

  write_message("Exported octree to " + path, 5);
}

const vector<int32_t>& OctreePartitioner::get_cells(
  Position r, Direction u) const
{
  if (!bounds_.contains(r)) {
    // discount this get_cells call from the stats to only focus on points
    // within the octree
    return get_cells_fallback(r, u);
  }

  // A normal octree traversal algorithm would store a node's center in the
  // nodes and use that to determine which child node to proceed to. However,
  // storing a node's center in a node is bad for the cache hit rate. As it
  // turns out, we can deduce the node's center on the fly by iteratively
  // changing the center by certain displacement every time we go down a level.

  // Initialize our node dimensions.
  Position node_dim;
  for (int i = 0; i < 3; i++) {
    node_dim[i] = (bounds_.max_[i] - bounds_.min_[i]) * 0.5;
  }

  Position center = bounds_.get_center();
  auto current = nodes_[0];
  while (!current.is_leaf()) {
    int idx = 0;
    for (int i = 0; i < 3; i++) {
      // Determine position relative to the center on the current axis.
      bool less = (r[i] < center[i]);
      // Update our node center
      node_dim[i] *= 0.5;
      center[i] += node_dim[i] * (less ? -1 : 1);
      // Update our index
      idx = 2 * idx + int(less);
    }

    // Proceed to next node
    current = nodes_[current.read_data() + idx];
  }

  const auto& cells = cell_data_[current.read_data()];

  if (cells.empty()) {
    // If we hit an empty cell, proceed to fallback directly
    return get_cells_fallback(r, u);
  }

  return cells;
}

const vector<int32_t>& OctreePartitioner::get_cells_fallback(
  Position r, Direction u) const
{
  return fallback.get_cells(r, u);
}

// This method is the exact same as the information refilling method in the
// octree constructor, the only difference being that it works on the compressed
// version of the octree instead the compressed one. Since it is pretty slow
// (and for some reason creates an octree that has the same failure rate but is
// slower), only use it if you are reading from a file (which would not have had
// any information refilling to cut down on file size).
void OctreePartitioner::refill_information()
{
  // Octree node but packed with extra info to make it easier to work with.
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

  // First, decompress the nodes to make them easier to work with
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

  // Forward propagation
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

  // Down propagation.
  std::stack<OctreeNodeExtraInfo*> uninfo_refilled_nodes;
  uninfo_refilled_nodes.push(nodes_to_propagate.back());
  while (!uninfo_refilled_nodes.empty()) {
    auto cur = uninfo_refilled_nodes.top();
    uninfo_refilled_nodes.pop();

    bool should_down_propagate = false;
    for (int i = 0; i < 8; i++) {
      if (cur->children_[i].is_leaf()) {
        should_down_propagate = true;
        break;
      }
    }

    auto collection_point = cur;
    for (int i = 0; i < INFORMATION_REFILLING_START_DEPTH_OFFSET; i++) {
      if (collection_point->parent_) {
        collection_point = collection_point->eiparent_;
      }
    }

    const auto& unique_cells = collection_point->cells_;

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