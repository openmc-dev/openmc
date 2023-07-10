#ifndef KDTREE_PARTITIONER_H
#define KDTREE_PARTITIONER_H

#include "universe.h"
#include "aabb.h"
#include "zplane_partitioner.h"
#include "octree_partitioner.h"
#include <vector>
#include <set>

namespace openmc {

    struct kdTreeNode {
        kdTreeNode* children;
        std::vector<int> cells;
    };

    class KdTree {
    public:

    private:
    
    };

};

#endif