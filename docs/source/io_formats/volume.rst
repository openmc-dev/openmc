.. _io_volume:

==================
Volume File Format
==================

**/**

:Attributes: - **samples** (*int*) -- Number of samples
             - **lower_left** (*double[3]*) -- Lower-left coordinates of
               bounding box
             - **upper_right** (*double[3]*) -- Upper-right coordinates of
               bounding box

**/cell_<id>/**

:Datasets: - **volume** (*double[2]*) -- Calculated volume and its uncertainty
             in cubic centimeters
           - **nuclides** (*char[][]*) -- Names of nuclides identified in the
             cell
           - **atoms** (*double[][2]*) -- Total number of atoms of each nuclide
             and its uncertainty
