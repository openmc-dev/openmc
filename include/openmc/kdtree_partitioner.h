#ifndef OPENMC_KDTREE_PARTITIONER_H
#define OPENMC_KDTREE_PARTITIONER_H

#include "octree_partitioner.h"
#include "partitioner_utils.h"
#include "universe.h"
#include "zplane_partitioner.h"
#include <vector>

namespace openmc {

struct KdTreeUncompressedNode {
  KdTreeUncompressedNode();
  bool is_leaf() const;

  AABB box;

  uint8_t split_axis;
  double split_location;
  KdTreeUncompressedNode* children;

  std::vector<int> cells;
};

struct KdTreeNode {
  uint32_t data;
  double split;

  bool is_leaf() const;
  uint32_t index() const;
};

struct KdTreeConstructionTask {
  KdTreeUncompressedNode* node;
  std::vector<CellPointUncompressed> points;
  uint32_t depth;

  KdTreeConstructionTask();
  KdTreeConstructionTask(KdTreeUncompressedNode* node,
    const std::vector<CellPointUncompressed>& points, uint32_t depth);
};

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
  AABB bounds;
  std::vector<KdTreeNode> nodes;
  std::vector<std::vector<int>> cell_data;
  ZPlanePartitioner fallback;
};

}; // namespace openmc

#endif // OPENMC_KDTREE_PARTITIONER_H