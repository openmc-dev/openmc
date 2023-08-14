// This header file automatically includes all partitioner types
// It also defines a DefaultPartitioner type
// Thus, if we want to change the default partitioner type to a new partitioner
// type, we just have to change it here instead of updating it everywhere
#ifndef OPENMC_PARTITIONERS_H
#define OPENMC_PARTITIONERS_H

#include "bin_grid_partitioner.h"
#include "kdtree_partitioner.h"
#include "octree_partitioner.h"
#include "universe.h"
#include "zplane_partitioner.h"

#define PARTITIONER_FALLBACK_ENABLED

namespace openmc {

using DefaultPartitioner = BinGridPartitioner;

enum PartitionerTypeID {
  Octree = 1,
  KdTree = 2,
  BinGrid = 3,
  ZPlane = 4,
};
}

#endif