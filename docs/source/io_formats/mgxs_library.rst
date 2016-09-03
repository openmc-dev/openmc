.. _io_mgxs_library:

========================================
Multi-Group Cross Section Library Format
========================================

OpenMC can be run in continuous-energy mode or multi-group mode, provided the
nuclear data is available.  In continuous-energy mode, the
``cross_sections.xml`` file contains necessary meta-data for each data set,
including the name and a file system location where the complete library
can be found.  In multi-group mode, the multi-group meta-data and the
nuclear data itself is contained within an ``mgxs.h5``.  This portion of
the manual describes the format of the multi-group data library required
to be used in the ``mgxs.h5`` file.

The multi-group library is provided in the HDF5_ format.  This library must
provide some meta-data about the library itself (such as the number of
groups and the group structure, etc.) as well as the actual cross section
data itself for each of the necessary nuclides or materials.

.. _HDF5: http://www.hdfgroup.org/HDF5/

.. _mgxs_lib_spec:

--------------------------
MGXS Library Specification
--------------------------

**/**

:Attributes: - **groups** (*int*) -- Number of energy groups
             - **group structure** (*double[]*) -- Monotonically increasing
               list of group boundaries, in units of MeV.  The length of this
               array should be the number of groups plus 1.
             - **inverse velocities** (*double[]*) -- Average inverse velocitiy
               for each of the groups in the library. This is an optional
               attribute.

**/<library name>/**

The data within <library name> contains the temperature-dependent multi-group
data for the nuclide or material that it represents.

:Attributes: - **awr** (*double*) -- The atomic weight ratio (optional, i.e. it
               is not meaningful for material-wise data)
             - **fissionable** (*int*) -- Whether the dataset is fissionable
               (1) or not (0).
             - **representation** (*char[]*) -- The method used to generate and
               represent the multi-group cross sections.  That is, whether they
               were generated with scalar flux weighting (or reduced to a
               similar representation) and thus are angle-independent, or if the
               data was generated with angular dependent fluxes and thus the
               data is angle-dependent.  Valid values are either "isotropic" or
               "angle".
             - **num-azimuthal** (*int*) -- Number of equal width angular bins
               that the azimuthal angular domain is subdivided if the
               `representation` attribute is "angle". This parameter is
               ignored otherwise.
            - **num-polar** (*int*) -- Number of equal width angular bins
               that the polar angular domain is subdivided if the
               `representation` attribute is "angle". This parameter is
               ignored otherwise.
            - **scatter-type** (*char[]*) -- The representation of the
              scattering angular distribution.  The options are either
              "legendre", "histogram", or "tabular".
            - **order** (*int*) -- Either the Legendre order, number of bins,
              or number of points (depending on the value of `scatter-type`)
              used to describe the angular distribution associated with each group-to-group transfer probability.

**/<library name>/kTs/**

:Datasets: - **<TTT>K** (*double*) -- kT values (in MeV) for each Temperature
             TTT (in Kelvin), rounded to the nearest integer

**/<library name>/<TTT>K/**

Temperature-dependent data, provided for temperature <TTT>K.

:Datasets: - **total** (*double[]* or *double[][][]*) -- Total cross section.
             This is a 1-D vector if `representation` is "isotropic", or a 3-D
             vector if `representation` is "angle" with dimensions of
             [groups, azimuthal, polar].
           - **absorption** (*double[]* or *double[][][]*) -- Absorption
             cross section.
             This is a 1-D vector if `representation` is "isotropic", or a 3-D
             vector if `representation` is "angle" with dimensions of
             [groups, azimuthal, polar].
           - **scatter matrix** (*double[][][]* or *double[][][][][]*) --
             Scattering moment matrices. This dataset requires the scattering
             moment matrices presented with the columns representing incoming
             group and rows representing the outgoing group. That is,
             down-scatter will be above the diagonal of the resultant matrix.
             This matrix is repeated for every Legendre order (in order of
             increasing orders) if `scatter-type` is "legendre"; otherwise, this
             matrix is repeated for every bin of the histogram or tabular
             representation.  Finally, if ``representation`` is "angle", the
             above is repeated for every azimuthal angle and every polar angle,
             in that order.
           - **multiplicity matrix** (*double[][]* or *double[][][][]*) --
             Scattering multiplicity matrices.
             This dataset provides the ratio of neutrons produced in scattering
             collisions to the neutrons which undergo scattering collisions;
             that is, the multiplicity provides the code with a scaling factor
             to account for neutrons being produced in (n,xn) reactions. This
             information is assumed isotropic and therefore is not repeated for
             every Legendre moment or histogram/tabular bin.  This matrix
             follows the same arrangement as described for the `scatter`
             dataset, with the exception of the data needed to provide the
             scattering type information.
             This dataset is optional, if it is not provided no multiplication
             (i.e., values of 1.0) will be assumed.
           - **fission** (*double[]* or *double[][][]*) -- Fission
             cross section.
             This is a 1-D vector if `representation` is "isotropic", or a 3-D
             vector if `representation` is "angle" with dimensions of
             [groups, azimuthal, polar].  This is only required if the data set
             is fissionable and fission-tallies are expected to be used.
           - **kappa-fission** (*double[]* or *double[][][]*) -- Kappa-Fission
             (energy-release from fission) cross section.
             This is a 1-D vector if `representation` is "isotropic", or a 3-D
             vector if `representation` is "angle" with dimensions of
             [groups, azimuthal, polar].  This is only required if the data set
             is fissionable and fission-tallies are expected to be used.
           - **chi** (*double[]* or *double[][][]*) -- Fission neutron energy
             spectra.
             This is a 1-D vector if `representation` is "isotropic", or a 3-D
             vector if `representation` is "angle" with dimensions of
             [groups, azimuthal, polar].  This is only required if the data set
             is fissionable and fission-tallies are expected to be used.
           - **nu-fission** (*double[]* to *double[][][][]*) -- Nu-Fission
             cross section.
             If **chi** is provided, then `nu-fission` has the same
             dimensionality as `fission`.  If **chi** is not provided, then
             the `nu-fission` data must represent the fission neutron energy
             spectra as well and thus will have one additional dimension
             for the outgoing energy group.  In this case, `nu-fission` has the
             same dimensionality as `multiplicity matrix`.