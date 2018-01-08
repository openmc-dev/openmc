.. _usersguide_cross_sections:

===========================
Cross Section Configuration
===========================

In order to run a simulation with OpenMC, you will need cross section data for
each nuclide or material in your problem. OpenMC can be run in continuous-energy
or multi-group mode.

In continuous-energy mode, OpenMC uses a native `HDF5
<https://support.hdfgroup.org/HDF5/>`_ format (see :ref:`io_nuclear_data`) to
store all nuclear data. If you have ACE format data that was produced with
NJOY_, such as that distributed with MCNP_ or Serpent_, it can be converted to
the HDF5 format using the :ref:`scripts_ace` script (or :ref:`using the Python
API <create_xs_library>`).  Several sources provide openly available ACE data as
described below and can be easily converted using the provided scripts. The
TALYS-based evaluated nuclear data library, TENDL_, is also available in ACE
format. In addition to tabulated cross sections in the HDF5 files, OpenMC relies
on :ref:`windowed multipole <windowed_multipole>` data to perform on-the-fly
Doppler broadening.

In multi-group mode, OpenMC utilizes an HDF5-based library format which can be
used to describe nuclide- or material-specific quantities.

---------------------
Environment Variables
---------------------

When :ref:`scripts_openmc` is run, it will look for several environment
variables that indicate where cross sections can be found. While the location of
cross sections can also be indicated through the :class:`openmc.Materials` class
(or in the :ref:`materials.xml <io_materials>` file), if you always use the same
set of cross section data, it is often easier to just set an environment
variable that will be picked up by default every time OpenMC is run. The
following environment variables are used:

:envvar:`OPENMC_CROSS_SECTIONS`
  Indicates the path to the :ref:`cross_sections.xml <io_cross_sections>`
  summary file that is used to locate HDF5 format cross section libraries if the
  user has not specified :attr:`Materials.cross_sections` (equivalently, the
  :ref:`cross_sections` in :ref:`materials.xml <io_materials>`).

:envvar:`OPENMC_MULTIPOLE_LIBRARY`
  Indicates the path to a directory containing windowed multipole data if the
  user has not specified :attr:`Materials.multipole_library` (equivalently, the
  :ref:`multipole_library` in :ref:`materials.xml <io_materials>`)

:envvar:`OPENMC_MG_CROSS_SECTIONS`
  Indicates the path to the an :ref:`HDF5 file <io_mgxs_library>` that contains
  multi-group cross sections if the user has not specified
  :attr:`Materials.cross_sections` (equivalently, the :ref:`cross_sections` in
  :ref:`materials.xml <io_materials>`).

To set these environment variables persistently, export them from your shell
profile (``.profile`` or ``.bashrc`` in bash_).

.. _bash: http://www.linuxfromscratch.org/blfs/view/6.3/postlfs/profile.html

--------------------------------
Continuous-Energy Cross Sections
--------------------------------

Using ENDF/B-VII.1 Cross Sections from NNDC
-------------------------------------------

The NNDC_ provides ACE data from the ENDF/B-VII.1 neutron and thermal scattering
sublibraries at room temperature processed using NJOY_. To use this data with
OpenMC, the :ref:`scripts_nndc` script can be used to automatically download and
extract the ACE data, fix any deficiencies, and create an HDF5 library:

.. code-block:: sh

    openmc-get-nndc-data

At this point, you should set the :envvar:`OPENMC_CROSS_SECTIONS` environment
variable to the absolute path of the file ``nndc_hdf5/cross_sections.xml``. This
cross section set is used by the test suite.

Using JEFF Cross Sections from OECD/NEA
---------------------------------------

The NEA_ provides processed ACE data from the JEFF_ library. To use this data
with OpenMC, the :ref:`scripts_jeff` script can be used to automatically
download and extract the ACE data, fix any deficiencies, and create an HDF5
library.

.. code-block:: sh

    openmc-get-jeff-data

At this point, you should set the :envvar:`OPENMC_CROSS_SECTIONS` environment
variable to the absolute path of the file ``jeff-3.2-hdf5/cross_sections.xml``.

Using Cross Sections from MCNP
------------------------------

OpenMC provides two scripts (:ref:`scripts_mcnp70` and :ref:`scripts_mcnp71`)
that will automatically convert ENDF/B-VII.0 and ENDF/B-VII.1 ACE data that is
provided with MCNP5 or MCNP6. To convert the ENDF/B-VII.0 ACE files
(``endf70[a-k]`` and ``endf70sab``) into the native HDF5 format, run the
following:

.. code-block:: sh

    openmc-convert-mcnp70-data /path/to/mcnpdata/

where ``/path/to/mcnpdata`` is the directory containing the ``endf70[a-k]``
files.

To convert the ENDF/B-VII.1 ACE files (the endf71x and ENDF71SaB libraries), use
the following script:

.. code-block:: sh

    openmc-convert-mcnp71-data /path/to/mcnpdata

where ``/path/to/mcnpdata`` is the directory containing the ``endf71x`` and
``ENDF71SaB`` directories.

.. _other_cross_sections:

Using Other Cross Sections
--------------------------

If you have a library of ACE format cross sections other than those listed above
that you need to convert to OpenMC's HDF5 format, the :ref:`scripts_ace` script
can be used. There are four different ways you can specify ACE libraries that
are to be converted:

1. List each ACE library as a positional argument. This is very useful in
   conjunction with the usual shell utilities (ls, find, etc.).
2. Use the ``--xml`` option to specify a pre-v0.9 cross_sections.xml file.
3. Use the ``--xsdir`` option to specify a MCNP xsdir file.
4. Use the ``--xsdata`` option to specify a Serpent xsdata file.

The script does not use any extra information from cross_sections.xml/ xsdir/
xsdata files to determine whether the nuclide is metastable. Instead, the
``--metastable`` argument can be used to specify whether the ZAID naming
convention follows the NNDC data convention (1000*Z + A + 300 + 100*m), or the
MCNP data convention (essentially the same as NNDC, except that the first
metastable state of Am242 is 95242 and the ground state is 95642).

.. _create_xs_library:

Manually Creating a Library from ACE files
------------------------------------------

.. currentmodule:: openmc.data

The scripts described above use the :mod:`openmc.data` module in the Python API
to convert ACE data and create a :ref:`cross_sections.xml <io_cross_sections>`
file. For those who prefer to use the API directly, the
:class:`openmc.data.IncidentNeutron` and :class:`openmc.data.ThermalScattering`
classes can be used to read ACE data and convert it to HDF5. For
continuous-energy incident neutron data, use the
:meth:`IncidentNeutron.from_ace` class method to read in an existing ACE file
and the :meth:`IncidentNeutron.export_to_hdf5` method to write the data to an
HDF5 file.

::

  u235 = openmc.data.IncidentNeutron.from_ace('92235.710nc')
  u235.export_to_hdf5('U235.h5')

If you have multiple ACE files for the same nuclide at different temperatures,
you can use the :meth:`IncidentNeutron.add_temperature_from_ace` method to
append cross sections to an existing :class:`IncidentNeutron` instance::

  u235 = openmc.data.IncidentNeutron.from_ace('92235.710nc')
  for suffix in [711, 712, 713, 714, 715, 716]:
      u235.add_temperature_from_ace('92235.{}nc'.format(suffix))
  u235.export_to_hdf5('U235.h5')

Similar methods exist for thermal scattering data:

::

  light_water = openmc.data.ThermalScattering.from_ace('lwtr.20t')
  for suffix in range(21, 28):
      light_water.add_temperature_from_ace('lwtr.{}t'.format(suffix))
  light_water.export_to_hdf5('lwtr.h5')

Once you have created corresponding HDF5 files for each of your ACE files, you
can create a library and export it to XML using the
:class:`openmc.data.DataLibrary` class::

  library = openmc.data.DataLibrary()
  library.register_file('U235.h5')
  library.register_file('lwtr.h5')
  ...
  library.export_to_xml()

At this point, you will have a ``cross_sections.xml`` file that you can use in
OpenMC.

.. hint:: The :class:`IncidentNeutron` class allows you to view/modify cross
          sections, secondary angle/energy distributions, probability tables,
          etc. For a more thorough overview of the capabilities of this class,
          see the :ref:`notebook_nuclear_data` example notebook.

Manually Creating a Library from ENDF files
-------------------------------------------

If you need to create a nuclear data library and you do not already have
suitable ACE files or you need to further customize the data (for example,
adding more temperatures), the :meth:`IncidentNeutron.from_njoy` and
:meth:`ThermalScattering.from_njoy` methods can be used to create data instances
by directly running NJOY. Both methods require that you pass the name of ENDF
file(s) that are passed on to NJOY. For example, to generate data for Zr-92::

  zr92 = openmc.data.IncidentNeutron.from_njoy('n-040_Zr_092.endf')

By default, data is produced at room temperature, 293.6 K. You can also specify
a list of temperatures that you want data at::

  zr92 = openmc.data.IncidentNeutron.from_njoy(
      'n-040_Zr_092.endf', temperatures=[300., 600., 1000.])

The :meth:`IncidentNeutron.from_njoy` method assumes you have an executable
named ``njoy`` available on your path. If you want to explicitly name the
executable, the ``njoy_exec`` optional argument can be used. Additionally, the
``stdout`` argument can be used to show the progress of the NJOY run.

To generate a thermal scattering file, you need to specify both an ENDF incident
neutron sub-library file as well as a thermal neutron scattering sub-library
file; for example::

  light_water = openmc.data.ThermalScattering.from_njoy(
      'neutrons/n-001_H_001.endf', 'thermal_scatt/tsl-HinH2O.endf')

Once you have instances of :class:`IncidentNeutron` and
:class:`ThermalScattering`, a library can be created by using the
``export_to_hdf5()`` methods and the :class:`DataLibrary` class as described in
:ref:`create_xs_library`.

Enabling Resonance Scattering Treatments
----------------------------------------

In order for OpenMC to correctly treat elastic scattering in heavy nuclides
where low-lying resonances might be present (see
:ref:`energy_dependent_xs_model`), the elastic scattering cross section at 0 K
must be present. To add the 0 K elastic scattering cross section to existing
:class:`IncidentNeutron` instance, you can use the
:meth:`IncidentNeutron.add_elastic_0K_from_endf` method which requires an ENDF
file for the nuclide you are modifying::

  u238 = openmc.data.IncidentNeutron.from_hdf5('U238.h5')
  u238.add_elastic_0K_from_endf('n-092_U_238.endf')
  u238.export_to_hdf5('U238_with_0K.h5')

With 0 K elastic scattering data present, you can turn on a resonance scattering
method using :attr:`Settings.resonance_scattering`.

.. note:: The process of reconstructing resonances and generating tabulated 0 K
          cross sections can be computationally expensive, especially for
          nuclides like U-238 where thousands of resonances are present. Thus,
          running the :meth:`IncidentNeutron.add_elastic_0K_from_endf` method
          may take several minutes to complete.

-----------------------
Windowed Multipole Data
-----------------------

OpenMC is capable of using windowed multipole data for on-the-fly Doppler
broadening. While such data is not yet available for all nuclides, an
experimental multipole library is available that contains data for 70
nuclides. To obtain this library, you can run :ref:`scripts_multipole` which
will download and extract it into a ``wmp`` directory. Once the library has been
downloaded, set the :envvar:`OPENMC_MULTIPOLE_LIBRARY` environment variable (or
the :attr:`Materials.multipole_library` attribute) to the ``wmp`` directory.

--------------------------
Multi-Group Cross Sections
--------------------------

Multi-group cross section libraries are generally tailored to the specific
calculation to be performed.  Therefore, at this point in time, OpenMC is not
distributed with any pre-existing multi-group cross section libraries.
However, if  obtained or generated their own library, the user
should set the :envvar:`OPENMC_MG_CROSS_SECTIONS` environment variable
to the absolute path of the file library expected to used most frequently.

For an example of how to create a multi-group library, see
:ref:`notebook_mg_mode_part_i`.

.. _NJOY: https://njoy.github.io/NJOY2016/
.. _NNDC: http://www.nndc.bnl.gov/endf/b7.1/acefiles.html
.. _NEA: http://www.oecd-nea.org
.. _JEFF: https://www.oecd-nea.org/dbforms/data/eva/evatapes/jeff_32/
.. _MCNP: http://mcnp.lanl.gov
.. _Serpent: http://montecarlo.vtt.fi
.. _TENDL: https://tendl.web.psi.ch/tendl_2015/tendl2015.html
