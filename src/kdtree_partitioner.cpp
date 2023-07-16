#include "openmc/kdtree_partitioner.h"
#include "openmc/error.h"
#include "openmc/partitioners.h"
#include <queue>
#include <stack>

namespace openmc {

KdTreeUncompressedNode::KdTreeUncompressedNode() : children_(nullptr) {}

// Check whether the current node is a leaf or not by checking if children_ are
// nullptr
bool KdTreeUncompressedNode::is_leaf() const
{
  return (children_ == nullptr);
}

// Compute the volume heuristic (VH)
double KdTreeUncompressedNode::compute_vh() const
{
  return box_.volume() * num_unique_cells_;
}

// If a cell is a leaf, then the last two bits will be on
constexpr uint32_t KD_TREE_LEAF_FLAG = (0b11 << 30);

// Checks whether the last two bits are on
bool KdTreeNode::is_leaf() const
{
  return ((data_ & KD_TREE_LEAF_FLAG) == KD_TREE_LEAF_FLAG);
}

// Extracts the first 30 bits out of data_ and returns them
uint32_t KdTreeNode::index() const
{
  uint32_t index = data_ & ~KD_TREE_LEAF_FLAG;
  return index;
}

// Gets the number of unique cells captured by the vector
// This method works by assuming the elements are sorted according to theirs
// cell IDs and then compares adjacent IDs
template<typename T>
uint32_t get_num_unique_cells(const std::vector<T>& v)
{
  int prev_cell = -1;
  uint32_t num_unique = 0;
  for (int i = 0; i < v.size(); i++) {
    int cell = v[i];
    if (cell != prev_cell) {
      num_unique++;
      prev_cell = cell;
    }
  }

  return num_unique;
}

// This work works exactly like get_num_unique_cells, expect this time it stores
// the cell list
template<typename T>
void build_cells_list(
  KdTreeUncompressedNode& node, const std::vector<T>& source)
{
  // Reserve memory
  node.cells_.reserve(node.num_unique_cells_);

  // Go through each point to find the unique cells
  int prev_cell = -1;
  for (int i = 0; i < source.size(); i++) {
    int cell = source[i];
    if (cell != prev_cell) {
      node.cells_.push_back(cell);
      prev_cell = cell;
    }
  }
}

double generate_split(
  KdTreeUncompressedNode& parent, KdTreeUncompressedNode children[2], int axis)
{
  // Calculate the midpoint, this will be where we split
  double midpoint = (parent.box_.max_[axis] + parent.box_.min_[axis]) * 0.5;

  // Create temporary children and assign them split boxes
  children[0].box_ = parent.box_;
  children[1].box_ = parent.box_;
  children[0].box_.max_[axis] = midpoint;
  children[1].box_.min_[axis] = midpoint;

  return midpoint;
}

// This function splits the box of a parent node in its midpoint on a given axis
// It then scores the children generated from that split, and if the children's
// VH score is better than the previous best, it replaces the previous best
// children with the new best children
void search_best_median_split(KdTreeUncompressedNode& parent,
  KdTreeUncompressedNode temp_children[2], int best_unique_cells[2],
  double& best_child_vh, int axis)
{
  // Clear any point memory from past lists
  for (int i = 0; i < 2; i++) {
    temp_children[i].cells_.clear();
  }

  // Create the split
  double midpoint = generate_split(parent, temp_children, axis);

  // Copy over points to each temporary child node
  for (const auto& p : parent.points_) {
    if (p.pos_[axis] < midpoint) {
      temp_children[0].cells_.push_back(p.cell_);
    } else {
      temp_children[1].cells_.push_back(p.cell_);
    }
  }

  // Calculate the children vh
  double cur_vh = 0.0;
  for (int i = 0; i < 2; i++) {
    temp_children[i].num_unique_cells_ =
      get_num_unique_cells(temp_children[i].cells_);
    cur_vh += temp_children[i].compute_vh();
  }

  // Update our best split so far if we find a split with a better vh
  if (cur_vh < best_child_vh) {
    best_child_vh = cur_vh;

    for (int i = 0; i < 2; i++) {
      best_unique_cells[i] = temp_children[i].num_unique_cells_;
    }

    // Store the split information
    // If we go through with the split, we will reuse this information when
    // actually generating the children nodes
    parent.split_axis_ = axis;
    parent.split_location_ = midpoint;
  }
}

KdTreePartitioner::KdTreePartitioner(
  const Universe& univ, const AABB& bounds, int32_t max_depth)
  : fallback_(univ), bounds_(bounds)
{
  write_message("Building kd-tree partitioner...", 5);

  Timer construction_timer;
  construction_timer.start();

  // Initialize our root
  KdTreeUncompressedNode root;
  root.box_ = bounds_;
  root.depth_ = 0;
  root.points_ =
    binned_point_search<CellPointUncompressed>(univ, fallback_, bounds_);

  // Pre-construction work: sort points by their cell IDs
  std::sort(root.points_.begin(), root.points_.end());
  root.num_unique_cells_ = get_num_unique_cells(root.points_);

  // Initialize some metadata
  uint32_t num_nodes = 1;
  uint32_t num_leaves = 0;

  // These temporary children are used as preallocated chunks of memory during
  // subdivision
  KdTreeUncompressedNode temp_children[2];
  for (int i = 0; i < 2; i++) {
    temp_children[i].cells_.reserve(root.points_.size());
  }

  // This will allocate nodes for us and manage the memory
  NodeAllocator<KdTreeUncompressedNode, 2> node_alloc;

  // Here is the subdivison loop
  // We recursively process nodes via a stack
  // For each node, we try to find the best split for it
  // If we find a split that results in a lower VH, we go with it
  // If we are not at max depth, we push the results of the split to the stack
  std::stack<KdTreeUncompressedNode*> unsubdivided_nodes;
  unsubdivided_nodes.push(&root);
  while (!unsubdivided_nodes.empty()) {
    auto& cur = *unsubdivided_nodes.top();
    unsubdivided_nodes.pop();

    // find best split
    double best_child_vh = DBL_MAX;
    int best_unique_cells[2];
    for (int i = 0; i < 3; i++) {
      search_best_median_split(
        cur, temp_children, best_unique_cells, best_child_vh, i);
    }

    // Subdivide if we found a better VH
    if (0.9 * best_child_vh < cur.compute_vh()) {
      num_nodes += 2;

      // Set up our children
      cur.children_ = node_alloc.allocate();
      generate_split(cur, cur.children_, cur.split_axis_);
      for (int i = 0; i < 2; i++) {
        cur.children_[i].num_unique_cells_ = best_unique_cells[i];
        cur.children_[i].depth_ = cur.depth_ + 1;
      }

      if (cur.depth_ + 1 == max_depth) {
        // We don't need to do anything too complex besides just building the
        // cells list
        for (int i = 0; i < 2; i++) {
          build_cells_list(cur.children_[i], temp_children[i].cells_);
          num_leaves += 2;
        }
      } else {
        // We need to copy over the points to their respective cell
        for (const auto& p : cur.points_) {
          if (p.pos_[cur.split_axis_] < cur.split_location_) {
            cur.children_[0].points_.push_back(p);
          } else {
            cur.children_[1].points_.push_back(p);
          }
        }

        // Now push to stack
        for (int i = 0; i < 2; i++) {
          unsubdivided_nodes.push(&cur.children_[i]);
        }
      }

    } else {
      // Since we didn't find a better VH, just turn the current node into a
      // leaf
      build_cells_list(cur, cur.points_);
      num_leaves++;
    }
  }

  // Now we begin the process of converting the uncompressed nodes to the
  // compressed nodes
  nodes_.reserve(num_nodes);
  cell_data_.reserve(num_leaves);

  // Process each node recursively, this time using a queue
  // We go through each node, and if it is a parent, we basically set up its
  // index to next child (we need a queue to predict this index ahead of time)
  // If it is a leaf, we store information in cell_data_
  std::queue<KdTreeUncompressedNode*> uncompressed_nodes;
  uncompressed_nodes.push(&root);
  while (!uncompressed_nodes.empty()) {
    auto& cur = *uncompressed_nodes.front();
    uncompressed_nodes.pop();

    // Compress
    KdTreeNode comp_node;
    if (cur.is_leaf()) {
      // Add information for the leaf
      comp_node.data_ = cell_data_.size();
      comp_node.data_ |= KD_TREE_LEAF_FLAG;

      cell_data_.push_back(std::move(cur.cells_));
    } else {
      // Add information for the parent
      comp_node.data_ = (cur.split_axis_ << 30);
      comp_node.data_ += nodes_.size() + 1 + uncompressed_nodes.size();
      comp_node.split_ = cur.split_location_;

      for (int i = 0; i < 2; i++) {
        uncompressed_nodes.push(&cur.children_[i]);
      }
    }

    // Now add to nodes array
    nodes_.push_back(comp_node);
  }
  write_message("Kd-tree construction completed in " +
                std::to_string(construction_timer.elapsed()) + " seconds.");
}

KdTreePartitioner::~KdTreePartitioner() {}

KdTreePartitioner::KdTreePartitioner(
  const Universe& univ, const AABB& bounds, hid_t file)
  : fallback_(univ), bounds_(bounds)
{
  // read the nodes
  hid_t nodes_group = open_group(file, "node_data");
  read_attr_int(nodes_group, "n_nodes", &num_nodes_);
  std::vector<int32_t> serialized_nodes;
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
    nodes_[i].data_ = serialized_nodes[2 * i + 0];
    nodes_[i].split_ = reinterpret_cast<float&>(serialized_nodes[2 * i + 1]);
    if (nodes_[i].is_leaf()) {
      num_leaves_++;

      std::vector<int> cells;
      cells.resize(nodes_[i].index());

      for (int& cell : cells) {
        cell = comp_cell_data.at(next_cell_index);
        next_cell_index++;
      }

      nodes_[i].data_ =
        KD_TREE_LEAF_FLAG | static_cast<uint32_t>(cell_data_.size());
      cell_data_.push_back(std::move(cells));
    }
  }
}

void KdTreePartitioner::export_to_hdf5(const std::string& path) const
{
  std::vector<int32_t> comp_cell_data;
  for (int i = 0; i < cell_data_.size(); i++) {
    for (int j = 0; j < cell_data_[i].size(); j++) {
      comp_cell_data.push_back(cell_data_[i][j]);
    }
  }

  hid_t file = file_open(path, 'w');

  // Write header
  write_attribute(file, "filetype", "partitioner");

  // write general partitioner information
  write_attribute(
    file, "part_type", static_cast<int>(PartitionerTypeID::KdTree));
  write_dataset(file, "bounds_max", bounds_.max_);
  write_dataset(file, "bounds_min", bounds_.min_);

  // write the nodes
  hid_t nodes_group = create_group(file, "node_data");
  write_attribute(nodes_group, "n_nodes", nodes_.size());
  // Serialize the data
  std::vector<int32_t> serialized_nodes(nodes_.size() * 2);
  for (int i = 0; i < nodes_.size(); i++) {
    auto node = nodes_[i]; // copy

    if (node.is_leaf()) {
      node.data_ = KD_TREE_LEAF_FLAG |
                   static_cast<uint32_t>(cell_data_[node.index()].size());
    }

    serialized_nodes[2 * i + 0] = reinterpret_cast<int32_t&>(node.data_);
    serialized_nodes[2 * i + 1] = reinterpret_cast<int32_t&>(node.split_);
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

//! Return the list of cells that could contain the given coordinates.
const std::vector<int32_t>& KdTreePartitioner::get_cells(
  Position r, Direction u) const
{
  // Immediately return fallback cells if outside bounds
  if (!bounds_.contains(r)) {
    return get_cells_fallback(r, u);
  }

  // Traverse through loop
  const KdTreeNode* cur = &nodes_[0];
  while (true) {
    uint32_t axis = (cur->data_ >> 30);
    if (axis == 3) {
      // This node is a leaf, exit
      break;
    } else {
      // Continue to the next node
      cur = &nodes_[cur->index() + uint32_t(r[axis] > cur->split_)];
    }
  }

  const auto& cells = cell_data_[cur->index()];

  // If we have an empty node, return fallback
  if (cells.empty()) {
    return get_cells_fallback(r, u);
  }

  return cells;
}

const std::vector<int32_t>& KdTreePartitioner::get_cells_fallback(
  Position r, Direction u) const
{
  return fallback_.get_cells(r, u);
}

}; // namespace openmc