.. _usersguide_voxel:

======================
Voxel Plot File Format
======================

**/filetype** (*char[]*)

    String indicating the type of file.

**/num_voxels** (*int[3]*)

    Number of voxels in the x-, y-, and z- directions.

**/voxel_width** (*double[3]*)

    Width of a voxel in centimeters.

**/lower_left** (*double[3]*)

    Cartesian coordinates of the lower-left corner of the plot.

**/data** (*int[][][]*)

    Data for each voxel that represents a material or cell ID.
