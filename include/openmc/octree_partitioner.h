#ifndef OCTREE_PARTITIONER_H
#define OCTREE_PARTITIONER_H

#include "universe.h"
#include "aabb.h"
#include "zplane_partitioner.h"
#include <vector>
#include <set>

namespace openmc {

using namespace std;

struct OctreeNode {
  OctreeNode();

  bool is_leaf() const;
  int get_containing_child_index(const vec3& pos) const;
  OctreeNode& get_containing_child(const vec3& pos) const;
  std::vector<AABB> subdivide(const AABB& parent);
  
  vec3 center;
  std::vector<int> cells;
  OctreeNode* children; // if this is a leaf, then children will be null
  int depth;
  int id;
};

class OctreeNodeAllocator {
public:
  OctreeNodeAllocator();
  ~OctreeNodeAllocator();
  OctreeNode* allocate();
private:
  std::vector<OctreeNode*> pools;
  int last_pool_next_index;
};

class OctreePartitioner : public UniversePartitioner {
public:
  explicit OctreePartitioner(const Universe& univ, int target_cells_per_node=6, int max_depth=6, const std::string& file_path="octree.bin");
  explicit OctreePartitioner(const Universe& univ, const std::string& file_path);
  virtual ~OctreePartitioner() override;

  //! Return the list of cells that could contain the given coordinates.
  virtual const vector<int32_t>& get_cells(Position r, Direction u) const override;
  virtual const vector<int32_t>& get_cells_fallback(Position r, Direction u) const override;

  void write_to_file(const std::string& file_path);
  void read_from_file(const std::string& file_path);
private:
  //! The root node is where the searching begins
  OctreeNode root;
  AABB bounds;

  int num_nodes;

  // fallback if octree doesn't work
  ZPlanePartitioner fallback;

  OctreeNodeAllocator node_alloc;
};

}

#endif