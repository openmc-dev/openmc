.. _devguide_voxel:

=====================================
Voxel Plot Binary File Specifications
=====================================

The current revision of the voxel plot binary file is 1.

**integer(4) n_voxels_x**

    Number of voxels in the x direction

**integer(4) n_voxels_y**

    Number of voxels in the y direction

**integer(4) n_voxels_z**

    Number of voxels in the z direction

**real(8) width_voxel_x**

    Width of voxels in the x direction

**real(8) width_voxel_y**

    Width of voxels in the y direction

**real(8) width_voxel_z**

    Width of voxels in the z direction

**real(8) lower_left_x**

    Lower left x point of the voxel grid

**real(8) lower_left_y**

    Lower left y point of the voxel grid

**real(8) lower_left_z**

    Lower left z point of the voxel grid

*do x = 1, n_voxels_x*
    *do y = 1, n_voxels_y*
        *do z = 1, n_voxels_z*

            **integer(4) id**

                Cell or material id number at this voxel center. Set to -1 when
                cell not_found.
