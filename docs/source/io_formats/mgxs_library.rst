.. _io_mgxs_library:

========================================
Multi-Group Cross Section Library Format
========================================

OpenMC can be run in continuous-energy mode or multi-group mode, provided the
nuclear data is available.  In continuous-energy mode, the
``cross_sections.xml`` file contains necessary meta-data for each dataset,
including the name and a file system location where the complete library
can be found.  In multi-group mode, the multi-group meta-data and the
nuclear data itself is contained within an ``mgxs.h5`` file.  This portion of
the manual describes the format of the multi-group data library required
to be used in the ``mgxs.h5`` file.

The multi-group library is provided in the HDF5_ format.  This library must
provide some meta-data about the library itself (such as the number of
energy groups, delayed groups, and the energy group structure, etc.) as
well as the actual cross section data itself for each of the necessary
nuclides or materials.

The current version of the multi-group library file format is 1.0.

.. _HDF5: http://www.hdfgroup.org/HDF5/

.. _mgxs_lib_spec:

--------------------------
MGXS Library Specification
--------------------------

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file;
               for this library it will be 'mgxs'.
             - **version** (*int[2]*) -- Major and minor version of the
               multi-group library file format.
             - **energy_groups** (*int*) -- Number of energy groups
             - **delayed_groups** (*int*) -- Number of delayed groups (optional)
             - **group structure** (*double[]*) -- Monotonically increasing
               list of group boundaries, in units of eV.  The length of this
               array should be the number of groups plus 1.

**/<library name>/**

The data within <library name> contains the temperature-dependent multi-group
data for the nuclide or material that it represents.

:Attributes: - **atomic_weight_ratio** (*double*) -- The atomic weight ratio
               (optional, i.e. it is not meaningful for material-wise data).
             - **fissionable** (*bool*) -- Whether the dataset is fissionable
               (True) or not (False).
             - **representation** (*char[]*) -- The method used to generate and
               represent the multi-group cross sections.  That is, whether they
               were generated with scalar flux weighting (or reduced to a
               similar representation) and thus are angle-independent, or if the
               data was generated with angular dependent fluxes and thus the
               data is angle-dependent.  Valid values are either "isotropic" or
               "angle".
             - **num_azimuthal** (*int*) -- Number of equal width angular bins
               that the azimuthal angular domain is subdivided if the
               `representation` attribute is "angle". This parameter is
               ignored otherwise.
             - **num_polar** (*int*) -- Number of equal width angular bins
               that the polar angular domain is subdivided if the
               `representation` attribute is "angle". This parameter is
               ignored otherwise.
             - **scatter_format** (*char[]*) -- The representation of the
               scattering angular distribution.  The options are either
               "legendre", "histogram", or "tabular".  If not provided, the
               default of "legendre" will be assumed.
             - **order** (*int*) -- Either the Legendre order, number of bins,
               or number of points (depending on the value of `scatter_format`)
               used to describe the angular distribution associated with each
               group-to-group transfer probability.
             - **scatter_shape** (*char[]*) -- The shape of the provided
               scatter and multiplicity matrix. The values provided are strings
               describing the ordering the scattering array is provided in
               row-major (i.e., C/C++ and Python) indexing. Valid values are
               "[Order][G][G']" or "[Order][G'][G]" where "G'" denotes the
               secondary/outgoing energy groups, "G" denotes the incoming
               energy groups, and "Order" is the angular distribution index.
               This value is not required; if not the default value of
               "[Order][G][G']" will be assumed.

**/<library name>/kTs/**

:Datasets:
           - **<TTT>K** (*double*) -- kT values (in eV) for each temperature
             TTT (in Kelvin), rounded to the nearest integer

**/<library name>/<TTT>K/**

Temperature-dependent data, provided for temperature <TTT>K.

:Datasets: - **total** (*double[]* or *double[][][]*) -- Total cross section.
             This is a 1-D vector if `representation` is "isotropic", or a 3-D
             vector if `representation` is "angle" with dimensions of
             [polar][azimuthal][groups].
           - **absorption** (*double[]* or *double[][][]*) -- Absorption
             cross section.
             This is a 1-D vector if `representation` is "isotropic", or a 3-D
             vector if `representation` is "angle" with dimensions of
             [groups][azimuthal][polar].
           - **fission** (*double[]* or *double[][][]*) -- Fission
             cross section.
             This is a 1-D vector if `representation` is "isotropic", or a 3-D
             vector if `representation` is "angle" with dimensions of
             [polar][azimuthal][groups].  This is only required if the dataset
             is fissionable and fission-tallies are expected to be used.
           - **kappa-fission** (*double[]* or *double[][][]*) -- Kappa-Fission
             (energy-release from fission) cross section.
             This is a 1-D vector if `representation` is "isotropic", or a 3-D
             vector if `representation` is "angle" with dimensions of
             [polar][azimuthal][groups].  This is only required if the dataset
             is fissionable and fission-tallies are expected to be used.
           - **chi** (*double[]* or *double[][][]*) -- Fission neutron energy
             spectra.
             This is a 1-D vector if `representation` is "isotropic", or a 3-D
             vector if `representation` is "angle" with dimensions of
             [polar][azimuthal][groups].  This is only required if the dataset
             is fissionable and fission-tallies are expected to be used.
           - **nu-fission** (*double[]* to *double[][][][]*) -- Nu-Fission
             cross section.
             If **chi** is provided, then `nu-fission` has the same
             dimensionality as `fission`.  If **chi** is not provided, then
             the `nu-fission` data must represent the fission neutron energy
             spectra as well and thus will have one additional dimension
             for the outgoing energy group.  In this case, `nu-fission` has the
             same dimensionality as `multiplicity matrix`.
           - **inverse-velocity** (*double[]* or *double[][][]*) --
             Average inverse velocity for each of the groups in the library.
             This dataset is optional. This is a 1-D vector if `representation`
             is "isotropic", or a 3-D vector if `representation` is "angle"
             with dimensions of [polar][azimuthal][groups].

**/<library name>/<TTT>K/scatter_data/**

Data specific to neutron scattering for the temperature <TTT>K

:Datasets: - **g_min** (*int[]* or *int[][][]*) --
             Minimum (most energetic) groups with non-zero values of
             the scattering matrix provided.  If `scatter_shape` is
             "[Order][G][G']" then `g_min` will describe the minimum values
             of "G'" for each "G"; if `scatter_shape` is "[Order][G'][G]"
             then `g_min` will describe the minimum values of "G" for each "G'".
             These group numbers use the standard
             ordering where the fastest neutron energy group is group 1 while
             the slowest neutron energy group is group G.
             The dimensionality of `g_min` is:
             `g_min[g]`, or `g_min[num_polar][num_azimuthal][g]`.
             The former is used when `representation` is "isotropic", and the
             latter when `representation` is "angle".
           - **g_max** (*int[]* or *int[][][]*) --
             Similar to `g_min`, except this dataset describes the maximum
             (least energetic) groups with non-zero values of
             the scattering matrix.
           - **scatter_matrix** (*double[]*) -- Flattened representation of the
             scattering moment matrices. The pre-flattened array corresponds to
             the shape provied in `scatter_shape`, but if `representation` is
             "angle" the dimensionality in `scatter_shape` is prepended by
             "[num_polar][num_azimuthal]" dimensions. The right-most energy
             group dimension will only include the entries between `g_min` and
             `g_max`.
             dimension has a dimensionality of `g_min` to `g_max`.
           - **multiplicity_matrix** (*double[]*) -- Flattened representation of
             the scattering moment matrices. This dataset provides the code with
             a scaling factor to account for neutrons being produced in (n,xn)
             reactions. This is assumed isotropic and therefore is not repeated
             for every Legendre moment or histogram/tabular bin. This dataset is
             optional, if it is not provided no multiplication (i.e., values of
             1.0) will be assumed.
             The pre-flattened array is shapes consistent with `scatter_matrix`
             except the "[Order]" dimension in `scatter_shape` is ignored since
             this data is assumed isotropic.
