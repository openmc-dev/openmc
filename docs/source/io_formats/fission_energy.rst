.. _usersguide_fission_energy:

==================================
Fission Energy Release File Format
==================================

This file is a compact HDF5 representation of the ENDF MT=1, MF=458 data (see
ENDF-102_ for details).  It gives the information needed to compute the energy
carried away from fission reactions by each reaction product (e.g. fragment
nuclei, neutrons) which depends on the incident neutron energy.  OpenMC is
distributed with one of these files under
data/fission_Q_data_endfb71.h5.  More files of this format can be created from
ENDF files with the
``openmc.data.write_compact_458_library`` function.  They can be read with the
``openmc.data.FissionEnergyRelease.from_compact_hdf5`` class method.

:Attributes: - **comment** (*char[]*) -- An optional text comment
             - **component order** (*char[][]*) -- An array of strings
               specifying the order each reaction product occurs in the data
               arrays.  The components use the 2-3 letter abbreviations
               specified in ENDF-102 e.g. EFR for fission fragments and ENP for
               prompt neutrons.

**/<nuclide name>/**
    Nuclides are named by concatenating their atomic symbol and mass number. For
    example, 'U235' or 'Pu239'.  Metastable nuclides are appended with an
    '_m' and their metastable number.  For example, 'Am242_m1'

:Datasets:
          - **data** (*double[][][]*) -- The energy release coefficients.  The
             first axis indexes the component type.  The second axis specifies
             values or uncertainties.  The third axis indexes the polynomial
             order.  If the data uses the Sher-Beck format, then the last axis
             will have a length of one and ENDF-102 should be consulted for
             energy dependence.  Otherwise, the data uses the Madland format
             which is a polynomial of incident energy.

             For example, if 'EFR' is given first in the **component order**
             attribute and the data uses the Madland format, then the energy
             released in the form of fission fragments at an incident energy
             :math:`E` is given by

             .. math::
                \text{data}[0, 0, 0] + \text{data}[0, 0, 1] \cdot E
                + \text{data}[0, 0, 2] \cdot E^2 + \ldots
              
             And its uncertainty is

             .. math::
                \text{data}[0, 1, 0] + \text{data}[0, 1, 1] \cdot E
                + \text{data}[0, 1, 2] \cdot E^2 + \ldots

.. _ENDF-102: http://www.nndc.bnl.gov/endfdocs/ENDF-102-2012.pdf
