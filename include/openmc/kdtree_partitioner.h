#ifndef OPENMC_KDTREE_PARTITIONER_H
#define OPENMC_KDTREE_PARTITIONER_H

#include "octree_partitioner.h"
#include "partitioner_utils.h"
#include "universe.h"
#include "zplane_partitioner.h"
#include <vector>

namespace openmc {

// This ndoe contains all the information you would want to know about a node
// It is not cache efficient to traverse and takes up extra memory that is not
// used during traversal
struct KdTreeUncompressedNode {
  KdTreeUncompressedNode();
  bool is_leaf() const;

  AABB box;

  uint8_t split_axis;
  double split_location;
  KdTreeUncompressedNode* children;

  std::vector<int> cells;
};

// This node contains the bare minimum information for traversal
// It is cut down from to help with better memory efficiency and improve the
// cache hit rate
struct KdTreeNode {
  uint32_t data;
  double split;

  bool is_leaf() const;
  uint32_t index() const;
};

// This serves the same purpose a std::tuple would but with a nicer format
struct KdTreeConstructionTask {
  KdTreeUncompressedNode* node;
  std::vector<CellPointUncompressed> points;
  uint32_t depth;

  KdTreeConstructionTask();
  KdTreeConstructionTask(KdTreeUncompressedNode* node,
    const std::vector<CellPointUncompressed>& points, uint32_t depth);
};

// The actual partitioner class
class KdTreePartitioner : public UniversePartitioner {
public:
  explicit KdTreePartitioner(const Universe& univ);
  virtual ~KdTreePartitioner() override;

  //! Return the list of cells that could contain the given coordinates.
  virtual const std::vector<int32_t>& get_cells(
    Position r, Direction u) const override;
  virtual const std::vector<int32_t>& get_cells_fallback(
    Position r, Direction u) const override;

private:
  // The bounds of the kd tree. Calling get_cell on points outside this region
  // will result in the fallback being used.
  AABB bounds;
  // This array contains all the nodes in the kd tree, with the root being at
  // index 0
  std::vector<KdTreeNode> nodes;
  // The cell lists for each leaf are stored here. The are stored seperately to
  // improve cache efficiency during traversal
  std::vector<std::vector<int>> cell_data;
  // Utilze the z-plane partitioner in a fallback, as it very rarely fails to
  // find the cell a particle is in.
  ZPlanePartitioner fallback;
};

}; // namespace openmc

#endif // OPENMC_KDTREE_PARTITIONER_H