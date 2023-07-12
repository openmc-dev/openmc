// This header file automatically includes all partitioner types 
// It also defines a DefaultPartitioner type
// Thus, if we want to change the default partitioner type to a new partitioner type, we just have to change it here instead of updating it everywhere
#ifndef OPENMC_PARTITIONERS_H
#define OPENMC_PARTITIONERS_H

#include "universe.h"
#include "zplane_partitioner.h"
#include "octree_partitioner.h"
#include "kdtree_partitioner.h"

#define PARTITIONER_FALLBACK_ENABLED

namespace openmc {

    using DefaultPartitioner = OctreePartitioner;

}

#endif