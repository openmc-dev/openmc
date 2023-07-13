#ifndef OCTREE_PARTITIONER_H
#define OCTREE_PARTITIONER_H

#include "universe.h"
#include "aabb.h"
#include "zplane_partitioner.h"
#include <vector>
#include <set>

namespace openmc {

  // this contains the bare minimum for traversal
  struct OctreeNode {
  public:
    OctreeNode();

    void mark_as_leaf();
    bool is_leaf() const;
    
    void store_data(uint32_t data);
    uint32_t read_data() const;

    int get_child_index(const vec3& r) const;
  private:
    uint32_t data;
  };

  struct OctreeUncompressedNode {
    OctreeUncompressedNode();

    bool is_leaf() const;
    void subdivide();

    bool contains(int cell) const;

    uint16_t depth;
    uint16_t num_original_cells;
    uint16_t num_unique_cells;

    AABB box;

    OctreeUncompressedNode* children; // if this is a leaf, then children will be null
    OctreeUncompressedNode* parent;

    std::vector<int32_t> cells;
  };

  class OctreePartitioner : public UniversePartitioner {
  public:
    explicit OctreePartitioner(const Universe& univ, int target_cells_per_node=6);
    explicit OctreePartitioner(const Universe& univ, const std::string& file_path);
    virtual ~OctreePartitioner() override;

    //! Return the list of cells that could contain the given coordinates.
    virtual const vector<int32_t>& get_cells(Position r, Direction u) const override;
    virtual const vector<int32_t>& get_cells_fallback(Position r, Direction u) const override;

    virtual void export_to_file(const std::string& file_path) const override;
  private:
    void compress(const OctreeUncompressedNode& root);
    void refill_information();

    AABB bounds;

    std::vector<OctreeNode> nodes;
    std::vector<std::vector<int32_t>> cell_data;

    // fallback if octree doesn't work
    ZPlanePartitioner fallback;

    uint32_t num_nodes;
    uint32_t num_leaves;
  };

}

#endif