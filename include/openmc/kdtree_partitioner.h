#ifndef KDTREE_PARTITIONER_H
#define KDTREE_PARTITIONER_H

#include "universe.h"
#include "aabb.h"
#include "zplane_partitioner.h"
#include "octree_partitioner.h"
#include "binning.h"
#include <vector>

namespace openmc {

    struct KdTreeUncompressedNode {
        KdTreeUncompressedNode();
        bool is_leaf() const;

        AABB box;

        KdTreeUncompressedNode* children;
        std::vector<int> cells;
    };

    struct KdTreeNode {
        AABB box;
        uint32_t data; 

        bool is_leaf() const;
        uint32_t index() const;
    };

    struct KdTreeConstructionTask {
        KdTreeUncompressedNode* node;
        std::vector<CellPointUncompressed> points;
        uint32_t depth;

        KdTreeConstructionTask();
        KdTreeConstructionTask(KdTreeUncompressedNode* node, const std::vector<CellPointUncompressed>& points, uint32_t depth);
    };

    class KdTreePartitioner : public UniversePartitioner {
    public:
        explicit KdTreePartitioner(const Universe& univ);
        virtual ~KdTreePartitioner() override;

        //! Return the list of cells that could contain the given coordinates.
        virtual const std::vector<int32_t>& get_cells(Position r, Direction u) const override;
        virtual const std::vector<int32_t>& get_cells_fallback(Position r, Direction u) const override;

    private:
        AABB bounds;
        std::vector<KdTreeNode> nodes;
        std::vector<std::vector<int>> cell_data;
        ZPlanePartitioner fallback;
    };

};

#endif