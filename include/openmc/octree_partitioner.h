#ifndef OCTREE_PARTITIONER_H
#define OCTREE_PARTITIONER_H

#include "universe.h"
#include "aabb.h"
#include "zplane_partitioner.h"
#include <vector>
#include <set>

namespace openmc {

using namespace std;

// this contains the bare minimum for traversal
struct OctreeNode {
  int32_t data;

  OctreeNode();
  bool is_leaf() const;
  
  int get_child_index(const vec3& r) const;

  void mark_leaf();
  int get_cells_index() const;
};

/*
struct OctreeNodeSuperCompressed {
  union {
    int32_t leaf_flag;

    struct {
      int32_t first_child_index;
    } parent_data;

    struct {
      uint16_t cell_offset;
      uint16_t num_cells;
    } child_data;
  };

  OctreeNode();
  bool is_leaf() const;

  void mark_leaf();
};
*/

struct OctreeUncompressedNode {
  OctreeUncompressedNode();

  bool is_leaf() const;
  int get_containing_child_index(const vec3& pos) const;
  OctreeUncompressedNode& get_containing_child(const vec3& pos) const;
  std::vector<AABB> subdivide();

  bool contains(int cell) const;

  uint32_t id;
  uint16_t depth;
  uint16_t num_original_cells;
  uint16_t num_unique_cells;

  AABB box;

  OctreeUncompressedNode* children; // if this is a leaf, then children will be null
  OctreeUncompressedNode* parent;

  std::vector<int> cells;
};


// unlike the compressed node, this contains extra information that might not be used in traversal
struct OctreeNodeSerialized {
  union {
    int32_t first_child_index; 
    int32_t num_contained_cells;      
  };

  void mark_leaf();
  bool is_leaf() const;

  int32_t get_num_cells() const;
};


class OctreePartitioner : public UniversePartitioner {
public:
  explicit OctreePartitioner(const Universe& univ, int target_cells_per_node=6, int max_depth=6, const std::string& file_path="octree.omcp");
  explicit OctreePartitioner(const Universe& univ, const std::string& file_path);
  virtual ~OctreePartitioner() override;

  //! Return the list of cells that could contain the given coordinates.
  virtual const vector<int32_t>& get_cells(Position r, Direction u) const override;
  virtual const vector<int32_t>& get_cells_fallback(Position r, Direction u) const override;

  void write_to_file(const std::string& file_path, const OctreeUncompressedNode& root) const;
  void read_from_file(const std::string& file_path);
  void compress(const OctreeUncompressedNode& root);
private:
  std::vector<OctreeNode> nodes;
  std::vector<std::vector<int>> cell_data;
  AABB bounds;

  std::vector<int> skip_traversal_grid;
  vec3 skip_traversal_dim;

  uint32_t num_nodes;

  // fallback if octree doesn't work
  ZPlanePartitioner fallback;

  const OctreeNode* get_containing_node(Position r,  const OctreeNode* start) const;
};

}

#endif