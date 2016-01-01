.. _usersguide_materials-out:

=========================
Materials-out File Format
=========================

If writing materials is requested, distributed compositions of all the materials
will be written to a HDF5 file named ````materials-out.h5````. 

**/mat-i/n_nuclides** (*int[]*)

    Number of nuclides in material *i*.

**/mat-i/n_instances** (*int[]*)

    Number of instances in material *i*.

**/mat-i/comps** (*double[]*)

    The compositions of all nuclides for all instances in material *i*.

