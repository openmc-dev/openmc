# OMCP File Format

The OpenMC Partitoner (OMCP) file format stores partitioners/acceleration structures for OpenMC. Objects such as octrees (and in the future, kd-trees and other partitioners) are stored in this file. All contents of a OMCP file are stored in a binary little endian format.

Every OMCP file opens up with 3 32-bit integers denoting the OMCP file format version. For instance, if the current version is 1.0.0, then the file would open up with the integers 1, 0, and 0. This is used to determine whether a certain OMCP file can be read by the current build of OpenMC.

The next value is another 32-bit integer, this time denoting what partitioner is stored in the file. The values are used for various partitioners (as of OMCP 1.0.0) are:
```c
enum OMCPStructureType {
    Octree = 1,
    KdTree = 2,               // not yet implemented in OMCP
    ZPlanePartitioner = 3,    // not yet implemented in OMCP
};
```

Based on the type of partitioner, the rest of the file format will be completely different, so I've broken down each partitioner into its own section.

## Octree

The octree consists of two parts:
* An array of serialized nodes
* An array of integers that represent the cells of the leaves of the octree (called the cell data here on out)

Here is a rough diagram of what the octree looks like in the file:

| Item | Bounds of the octree |Number of nodes | Serialized nodes array | Size of cell data array (4 bytes) | Cell data array |
| --- | --- | --- | --- | --- | --- |
| Size (in bytes) | 24 | 4 | number of nodes * sizeof(OctreeNodeSerialized) | 4 | size of cell data array * 2 |

The bounds have the following layout on disk:
```cpp
// include/openmc/AABB.h
struct AABB {
    Position min;
	Position max;
};

```

A serialized node has the following layout on disk:
```cpp
// include/openmc/OctreePartitioner.h
struct OctreeNodeSerialized {
  union {
    int32_t first_child_index; 
    int32_t num_contained_cells;      
  };
};
```
The node at index 0 of the serialized nodes array is the root. 

To determine whether a node is a leaf or parent, the 32nd bit will be set to 1 if the node is a leaf. Otherwise, it is a parent. 

If a node is a parent, then `first_child_index` will be an index that points to a location in the serialized nodes array where its first child is. All children of a node are contiguous in memory. The ordering of the children nodes will be discussed later. 

If a node is a leaf, then `num_contained_cells` will store the number of cells in that leaf. All of a leaf's cells are stored contiguously in the cell data array. Notably, `OctreeNodeSerialized` does not give any indication where in the cell data array a particular leaf's cells start. To reduce file size, the cells of the leaves with the smallest node IDs come first (the ID of any particular node is equal to its index). Thus, you can implictly determine where a particular leaf's cells are.

### Children node ordering

Since there are 8 children nodes, the indicies of all nodes will be in {0b000, ..., 0b111}. If we set the value of each bit to 1 if the center of a child node's box is less than the center of parent node's box on each axis and 0 otherwise, we have a way to order the nodes. In OMCP, the axis ordering is XYZ (e.g. if X is less but Y and Z are more, then our index will be 0b100 in binary, or 4 in base 10). 

Credit goes to Gavin for this idea.

### Node Boxes

OMCP does not store the box of each node. Instead, the box of a node can be derived via the bounds. The box of the root can be set to the bounds, and then the box can be subdivided. Then, the subdivided boxes can be assigned to the root's children, and this can occur recursively for each parent node. 

## Kd Trees

Kd trees are not implemented yet in OpenMC, so in OMCP 1.0.0, there is no specification for reading/writing Kd trees. The value of `OMCPStructureType::KdTree` is reserved for future use.

## Z-Plane Partitioner

Although z-plane partitioners are implemented in OpenMC, OMCP 1.0.0 does not have a specification for reading/writing them. The value of `OMCPStructureType::ZPlanePartitioner` is reserved for future use.