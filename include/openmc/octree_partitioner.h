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
  vec3 center;
  int data;

  // a really hacking compression/better cache usage approach would be to do away with vectors and do something like this
  // union {
  //    vec3 center;
  //    struct {
  //       int* cells;
  //       int num_cells;
  //    };
  // };
  // then we wouldn't have to follow get_cells_index to get our cells, we could just use the pointer for better cache

  OctreeNode();
  bool is_leaf() const;
  
  int get_child_index(const vec3& r) const;

  void mark_leaf();
  int get_cells_index() const;
};

struct OctreeUncompressedNode {
  OctreeUncompressedNode();

  bool is_leaf() const;
  int get_containing_child_index(const vec3& pos) const;
  OctreeUncompressedNode& get_containing_child(const vec3& pos) const;
  std::vector<AABB> subdivide(const AABB& parent);

  bool contains(int cell) const;
  
  vec3 center;
  AABB box;
  std::vector<int> cells;
  OctreeUncompressedNode* children; // if this is a leaf, then children will be null
  int depth;
  int id;
};


// unlike the compressed node, this contains extra information that might not be used in traversal
struct OctreeNodeSerialized {
    int id;

    vec3 center;
    bool is_leaf;

    // parent data
    int first_child_index;

    // leaf data
    int contained_cells_index;
    int num_contained_cells;
};


class OctreePartitioner : public UniversePartitioner {
public:
  explicit OctreePartitioner(const Universe& univ, int target_cells_per_node=6, int max_depth=6, const std::string& file_path="octree.bin");
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

  int num_nodes;

  // fallback if octree doesn't work
  ZPlanePartitioner fallback;
};

}

#endif