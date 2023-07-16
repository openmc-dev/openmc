#ifndef OPENMC_OCTREE_PARTITIONER_H
#define OPENMC_OCTREE_PARTITIONER_H

#include "partitioner_utils.h"
#include "universe.h"
#include "zplane_partitioner.h"
#include <set>
#include <vector>

namespace openmc {

// This struct contains the bare minimum information for traversal
// A leaf flag is stored in the 32nd bit
// The remaining 31 bits are dedicated to parent or node specfic information
// In the case of a parent, it stores the index of the first child node in the
// node array In the case of a leaf, it stores the index of the cell list in the
// cell_list_ array
struct OctreeNode {
public:
  OctreeNode();

  // Set 32nd bit to on
  void mark_as_leaf();
  // Check if 32nd bit is on
  bool is_leaf() const;

  // Store a 31 bit integer in the first 31 bits of data_
  void store_data(uint32_t data);
  // Extract the first 31 bits of data_
  uint32_t read_data() const;

private:
  uint32_t data_;
};

// This struct contains extra information about a node that is not needed for
// traversal
struct OctreeUncompressedNode {
  OctreeUncompressedNode();

  // Checks whether children_ is nullptr
  bool is_leaf() const;
  // Initialize children nodes
  void subdivide();
  // Check whether cell is in cell_ (assumes cells_ is sorted)
  bool contains(int cell) const;

  // Bounding box of node
  AABB box_;
  // List of cells contained within this node
  std::vector<int32_t> cells_;

  // Pointer to first child node (all children nodes are stored contigously in
  // memory). If the node is a leaf, then children will be null,
  OctreeUncompressedNode* children_;

  // Pointer to parent node. Used in information refilling.
  OctreeUncompressedNode* parent_;

  // Depth of current node
  uint16_t depth_;
  // Number of cells contained (equal to cells_.size() but unlike cells_.size(),
  // you don't have to build the cells list to get the number of unique cells)
  // Instead you can store them in an integer. This is useful for saving
  // performance.
  uint16_t num_unique_cells_;
  // Number of cells contained in cells_ before information refilling
  uint16_t num_original_cells_;
};

class OctreePartitioner : public UniversePartitioner {
public:
  explicit OctreePartitioner(
    const Universe& univ, const AABB& bounds, int target_cells_per_node = 6);
  explicit OctreePartitioner(
    const Universe& univ, const std::string& file_path);
  virtual ~OctreePartitioner() override;

  //! Return the list of cells that could contain the given coordinates.
  virtual const vector<int32_t>& get_cells(
    Position r, Direction u) const override;
  virtual const vector<int32_t>& get_cells_fallback(
    Position r, Direction u) const override;

  virtual void export_to_file(const std::string& file_path) const override;

private:
  void compress(const OctreeUncompressedNode& root);
  void refill_information();

  AABB bounds_;

  std::vector<OctreeNode> nodes_;
  std::vector<std::vector<int32_t>> cell_data_;

  // fallback if octree doesn't work
  ZPlanePartitioner fallback;

  // meta data, only used when exporting to file
  uint32_t num_nodes_;
  uint32_t num_leaves_;
  std::vector<uint32_t> orig_size_;
};

} // namespace openmc

#endif // OPENMC_OCTREE_PARTITIONER_H