.. _io_formats_partitioners.rst:

=======================
Partitioner file format
=======================

The partitioner file format will only have a subset of the following attributes. The exact subset is determined by the partitioner type.

**/**

:Datasets: - **bounds_min** (*Position*) -- the minimum bounds of the partitioner.
           - **bounds_max** (*Position*) -- the maximum bounds of the partitioner.

:Attributes: - **filetype** (*char[]*) -- A file type containing the string "partitioner".
             - **part_type** (*enum PartitionerTypeID*) -- an enum to corresponds to the partitioner type.

**/bin_grid/**

:Attributes: - **grid_res** (*int*) -- How many bins are on each axis.

**/bin_grid/bins/**

:Datasets: - **{i}** (*int[]*) -- The cells within the *i*th bin, where *i* is in {0...grid_res^3 - 1}

**/kdtree/node_data/**

:Datasets: - **nodes** (*int[]*) -- An array containing serialized nodes. Each node is a two-int pair, the first containing node-specific data and the second being the split location reinterpreted as a float.

:Attributes: - **n_nodes** (*size_t*) -- The number of nodes in the KD tree.

**/kdtree/cell_data/**

:Datasets: - **cells** (*int[]*) -- The list of cells of all child nodes combined into a single array. Leaf nodes with smaller IDs always have their cell lists come earlier in the array.

:Attributes: - **n_cells** (*size_t*) -- The number of cells in the KD tree.

**/octree/node_data/**

:Datasets: - **nodes** (*int[]*) -- An array containing serialized nodes.

:Attributes: - **n_nodes** (*size_t*) -- The number of nodes in the KD tree.

**/octree/cell_data/**

:Datasets: - **cells** (*int[]*) -- The list of cells of all child nodes combined into a single array. Leaf nodes with smaller IDs always have their cell lists come earlier in the array.

:Attributes: - **n_cells** (*size_t*) -- The number of cells in the KD tree.

**/zplane/**

:Attributes: - **surfs** (*int[]*) -- The list of surface IDs that partition the universe.

**/zplane/partitions/**

:Datasets: - **partition{i}** (*int[]*) -- The list of cells contained within partition *i*, where *i* is in {0...n_partitions - 1}.

:Attributes: - **n_partitions** (*int*) -- The number of partitions.


