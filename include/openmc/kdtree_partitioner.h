#ifndef OPENMC_KDTREE_PARTITIONER_H
#define OPENMC_KDTREE_PARTITIONER_H

#include <vector>

#include "octree_partitioner.h"
#include "partitioner_utils.h"
#include "universe.h"
#include "zplane_partitioner.h"

namespace openmc {

// This ndoe contains all the information you would want to know about a node
// It is not cache efficient to traverse and takes up extra memory that is not
// used during traversal
struct KdTreeUncompressedNode {
  KdTreeUncompressedNode();

  // This returns where or not children_ is nullptr
  bool is_leaf() const;

  // Compute the volume heuristic (VH)
  double compute_vh() const;

  // Bounding box of the node
  AABB box_;

  // If this node is a parent, then split_axis_ marks which axis (x=0, y=1, z=2)
  // the split is on
  uint8_t split_axis_;
  // Split location marks what point along the split axis the node is split on
  double split_location_;

  // Location of the first child if the current node is a parent (for all
  // tree-based partitioners, nodes are stored contigously in memory) If the
  // current node is a leaf, this will be nullptr
  KdTreeUncompressedNode* children_;

  // The list of unique cells for this node.
  std::vector<int> cells_;

  // The number of unique cells
  // Storing this allows us to know how many unique cells are in points_ without
  // having to build cells_
  uint32_t num_unique_cells_;

  // This contains the randomly sampled points contained within this node
  std::vector<CellPointUncompressed> points_;

  // The depth of the node in the tree
  uint32_t depth_;
};

// This node contains the bare minimum information for traversal
// It is cut down from to help with better memory efficiency and improve the
// cache hit rate
struct KdTreeNode {
  // Last two bits pack together the leaf flag and what axis the node was split
  // on if it is a parent If the last two bits equal 0b00, 0b01, or 0b10, then
  // it is a parent an the value of the bits correspond to the split axis If the
  // the bits equal 0b11, then it is a leaf The remaining 30 bits correspond to
  // what index in KdTreePartitioner::cell_data_ is the cell list if the node is
  // a leaf If the node is a parent, then the 30 bits correspond to the index in
  // KdTreePartitioner::nodes where the first child is
  uint32_t data_;
  // If this node is a parent, then this marks where along the axis the node is
  // split I use single precision here in the interest of cache efficiency while
  // traversal
  float split_;

  // Returns whether the current node is a leaf or not
  bool is_leaf() const;

  // This returns the data_ but with the flags/split axis information bits set
  // to zero In the case of a parent, this will return the index to the first
  // child node In the case of a leaf, this will return an index to cell_data_
  // to mark where its cell list is located
  uint32_t index() const;
};

// The actual partitioner class
class KdTreePartitioner : public UniversePartitioner {
public:
  explicit KdTreePartitioner(const Universe& univ, uint32_t max_depth);
  virtual ~KdTreePartitioner() override;

  //! Return the list of cells that could contain the given coordinates.
  virtual const std::vector<int32_t>& get_cells(
    Position r, Direction u) const override;
  virtual const std::vector<int32_t>& get_cells_fallback(
    Position r, Direction u) const override;

private:
  // The bounds of the kd tree. Calling get_cell on points outside this region
  // will result in the fallback being used.
  AABB bounds_;
  // This array contains all the nodes in the kd tree, with the root being at
  // index 0
  std::vector<KdTreeNode> nodes_;
  // The cell lists for each leaf are stored here. The are stored seperately to
  // improve cache efficiency during traversal
  std::vector<std::vector<int>> cell_data_;
  // Utilze the z-plane partitioner in a fallback, as it very rarely fails to
  // find the cell a particle is in.
  ZPlanePartitioner fallback_;
};

}; // namespace openmc

#endif // OPENMC_KDTREE_PARTITIONER_H