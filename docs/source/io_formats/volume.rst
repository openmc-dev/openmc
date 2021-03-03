.. _io_volume:

==================
Volume File Format
==================

The current version of the volume file format is 1.0.

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file.
             - **version** (*int[2]*) -- Major and minor version of the summary
               file format.
             - **openmc_version** (*int[3]*) -- Major, minor, and release
               version number for OpenMC.
             - **git_sha1** (*char[40]*) -- Git commit SHA-1 hash.
             - **date_and_time** (*char[]*) -- Date and time the summary was
               written.
             - **domain_type** (*char[]*) -- The type of domain for which
               volumes are calculated, either 'cell', 'material', or 'universe'.
             - **samples** (*int*) -- Number of samples
             - **lower_left** (*double[3]*) -- Lower-left coordinates of
               bounding box
             - **upper_right** (*double[3]*) -- Upper-right coordinates of
               bounding box
             - **threshold** (*double*) -- Threshold used for volume uncertainty
             - **trigger_type** (*char[]*) -- Trigger type used for volume uncertainty

**/domain_<id>/**

:Datasets: - **volume** (*double[2]*) -- Calculated volume and its uncertainty
             in cubic centimeters
           - **nuclides** (*char[][]*) -- Names of nuclides identified in the
             domain
           - **atoms** (*double[][2]*) -- Total number of atoms of each nuclide
             and its uncertainty
