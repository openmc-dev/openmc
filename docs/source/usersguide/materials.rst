.. _usersguide_materials:

.. currentmodule:: openmc

=====================
Material Compositions
=====================

Materials in OpenMC are defined as a set of nuclides/elements at specified
densities and are created using the :class:`openmc.Material` class. Once a
material has been instantiated, nuclides can be added with
:meth:`Material.add_nuclide` and elements can be added with
:meth:`Material.add_element`. Densities can be specified using atom fractions or
weight fractions. For example, to create a material and add Gd152 at 0.5 atom
percent, you'd run::

   mat = openmc.Material()
   mat.add_nuclide('Gd152', 0.5, 'ao')

The third argument to :meth:`Material.add_nuclide` can also be 'wo' for weight
percent. The densities specified for each nuclide/element are relative and are
renormalized based on the total density of the material. The total density is
set using the :meth:`Material.set_density` method. The density can be specified
in gram per cubic centimeter ('g/cm3'), atom per barn-cm ('atom/b-cm'), or
kilogram per cubic meter ('kg/m3'), e.g.,

::

   mat.set_density('g/cm3', 4.5)

----------------
Natural Elements
----------------

The :meth:`Material.add_element` method works exactly the same as
:meth:`Material.add_nuclide`, except that instead of specifying a single isotope
of an element, you specify the element itself. For example,

::

   mat.add_element('C', 1.0)

This method can also accept case-insensitive element names such as

::

  mat.add_element('aluminium', 1.0)

Internally, OpenMC stores data on the atomic masses and natural abundances of
all known isotopes and then uses this data to determine what isotopes should be
added to the material. When the material is later exported to XML for use by the
:ref:`scripts_openmc` executable, you'll see that any natural elements were
expanded to the naturally-occurring isotopes.

The :meth:`Material.add_element` method can also be used to add uranium at a
specified enrichment through the `enrichment` argument. For example, the
following would add 3.2% enriched uranium to a material::

   mat.add_element('U', 1.0, enrichment=3.2)

In addition to U235 and U238, concentrations of U234 and U236 will be present
and are determined through a correlation based on measured data.

It is also possible to perform enrichment of any element that is composed
of two naturally-occurring isotopes (e.g., Li or B) in terms of atomic percent.
To invoke this, provide the additional argument `enrichment_target` to
:meth:`Material.add_element`. For example the following would enrich B10
to 30ao%::

   mat.add_element('B', 1.0, enrichment=30.0, enrichment_target='B10')

In order to enrich an isotope in terms of mass percent (wo%), provide the extra
argument `enrichment_type`. For example the following would enrich Li6 to 15wo%::

   mat.add_element('Li', 1.0, enrichment=15.0, enrichment_target='Li6',
                   enrichment_type='wo')

Some older cross section libraries do not have all naturally-occurring isotopes
for a given element. For example, in ENDF/B-VII.1, cross section evaluations are
given for O16 and O17 but not for O18. If OpenMC is aware of what cross sections
you will be using through ``openmc.config['cross_sections']``, it will attempt
to only put isotopes in your model for which you have cross section data. In the
case of oxygen in ENDF/B-VII.1, the abundance of O18 would end up being lumped
with O16.

------------------
Material Libraries
------------------

OpenMC includes material compositions from Revision 2 of the `PNNL Compendium of
Material Composition Data for Radiation Transport Modeling
<https://doi.org/10.2172/1782721>`_. A material can be created from this library
with :meth:`Material.from_library`, for example::

   sodium_oxide = openmc.Material.from_library('Sodium Oxide')

Keyword arguments accepted by :class:`Material` can also be specified, such as
a material ID or temperature::

   salt_water = openmc.Material.from_library(
       'Salt Water', material_id=10, temperature=293.15)

Material names are case sensitive, but whitespace is normalized when looking up
a name. The default bundled library is ``'pnnl_v2'``. Each call returns a new,
independent material with the density and composition reported in the compendium.
The available material names can be retrieved without creating any materials::

   names = openmc.Material.get_library_material_names()

Natural compositions in the PNNL library are stored as elemental atom fractions
and expanded using :meth:`Material.add_element`. As described above, OpenMC will
account for the nuclides available in the cross section library indicated by
``openmc.config['cross_sections']``. Materials that are explicitly isotopic in
the compendium, such as enriched uranium or He-3 proportional gas, retain those
specified nuclides. OpenMC cannot substitute for a required isotope if it is
absent from the selected nuclear data library.

The compendium values are representative material definitions; actual density
and composition can vary. Users should confirm that a library material is
appropriate for their application.

Custom libraries can be registered from local JSON files. Registration applies
to the current Python process and cannot replace a library with the same name::

   openmc.Material.register_library('custom', 'custom_materials.json')
   names = openmc.Material.get_library_material_names('custom')
   material = openmc.Material.from_library(
       'My Material', library='custom', temperature=600.0)

Custom files use the same schema as the bundled library. A minimal library
containing an elemental composition is shown below:

.. code-block:: json

   {
     "schema_version": 1,
     "density_units": "g/cm3",
     "percent_type": "ao",
     "materials": {
       "My Material": {
         "density": 1.0,
         "elements": {
           "H": 0.666667,
           "O": 0.333333
         }
       }
     }
   }

Each material can contain ``elements``, ``nuclides``, or both. Elemental
components are expanded using :meth:`Material.add_element` and
``openmc.config['cross_sections']``. Nuclide components retain their explicitly
specified isotopic composition. The top-level ``percent_type`` and
``density_units`` values are passed to :class:`Material`, which validates them
along with the density and components when a material is created. Component
values can be relative amounts and do not need to sum to one. Additional
provenance metadata, such as a top-level ``source`` or a per-material note, is
permitted and ignored when constructing a material.

-----------------------
Thermal Scattering Data
-----------------------

If you have a moderating material in your model like water or graphite, you
should assign thermal scattering data (so-called :math:`S(\alpha,\beta)`) using
the :meth:`Material.add_s_alpha_beta` method. For example, to model light water,
you would need to add hydrogen and oxygen to a material and then assign the
``c_H_in_H2O`` thermal scattering data::

   water = openmc.Material()
   water.add_nuclide('H1', 2.0)
   water.add_nuclide('O16', 1.0)
   water.add_s_alpha_beta('c_H_in_H2O')
   water.set_density('g/cm3', 1.0)

.. _usersguide_naming:

-------------------------
Adding NCrystal materials
-------------------------

Additional support for thermal scattering can be added by using NCrystal_. The
:meth:`Material.from_ncrystal` class method generates a :class:`openmc.Material`
object from an `NCrystal configuration string
<https://github.com/mctools/ncrystal/wiki/Using-NCrystal#uniform-material-configuration-syntax>`_.
Temperature, material composition, and density are passed from the configuration
string and the `NCMAT file
<https://github.com/mctools/ncrystal/wiki/NCMAT-format>`_ that define the
material, e.g.::

   mat = openmc.Material.from_ncrystal('Al_sg225.ncmat;temp=300K')

defines a material containing polycrystalline alumnium,

::

   mat = openmc.Material.from_ncrystal("""Ge_sg227.ncmat;dcutoff=0.5;mos=40arcsec;
                                           dir1=@crys_hkl:5,1,1@lab:0,0,1;
                                           dir2=@crys_hkl:0,-1,1@lab:0,1,0""")

defines an oriented germanium single crystal with 40 arcsec mosaicity.

NCrystal only handles low energy neutron interactions. Other interactions are
provided by standard ACE files. NCrystal_ comes with a `predefined library
<https://github.com/mctools/ncrystal/wiki/Data-library>`_ but more materials can
be added by creating NCMAT files or on-the-fly in the configuration string.

.. warning:: Currently, NCrystal_ materials cannot be modified after they are created.
             Density, temperature and composition should be defined in the
             configuration string or the NCMAT file.

.. _NCrystal: https://github.com/mctools/ncrystal

------------------
Naming Conventions
------------------

OpenMC uses the GNDS_ naming convention for nuclides, metastable states, and
compounds:

:Nuclides: ``SymA`` where "A" is the mass number (e.g., ``Fe56``)
:Elements: ``Sym0`` (e.g., ``Fe0`` or ``C0``)
:Excited states: ``SymA_eN`` (e.g., ``V51_e1`` for the first excited state of
                 Vanadium-51.) This is only used in decay data.
:Metastable states: ``SymA_mN`` (e.g., ``Am242_m1`` for the first excited state
                    of Americium-242).
:Compounds: ``c_String_Describing_Material`` (e.g., ``c_H_in_H2O``). Used for
            thermal scattering data.

.. important:: The element syntax, e.g., ``C0``, is only used when the cross
               section evaluation is an elemental evaluation, like carbon in
               ENDF/B-VII.1! If you are adding an element via
               :meth:`Material.add_element`, just use ``Sym``.

.. _GNDS: https://www.oecd-nea.org/jcms/pl_39689/specifications-for-the-generalised-nuclear-database-structure-gnds

-----------
Temperature
-----------

Some Monte Carlo codes define temperature implicitly through the cross section
data, which is itself given only at a particular temperature. In OpenMC, the
material definition is decoupled from the specification of temperature. Instead,
temperatures are assigned to :ref:`cells <usersguide_cells>`
directly. Alternatively, a default temperature can be assigned to a material
that is to be applied to any cell where the material is used. In the absence of
any cell or material temperature specification, a global default temperature can
be set that is applied to all cells and materials. Anytime a material
temperature is specified, it will override the global default
temperature. Similarly, anytime a cell temperatures is specified, it will
override the material or global default temperature. All temperatures should be
given in units of Kelvin.

To assign a default material temperature, one should use the ``temperature``
attribute, e.g.,

::

   hot_fuel = openmc.Material()
   hot_fuel.temperature = 1200.0  # temperature in Kelvin

.. warning:: MCNP_ users should be aware that OpenMC does not use the concept of
             cross section suffixes like "71c" or "80c". Temperatures in Kelvin
             should be assigned directly per material or per cell using the
             :attr:`Material.temperature` or :attr:`Cell.temperature`
             attributes, respectively.

-----------------
Material Mixtures
-----------------

In OpenMC it is possible to mix any number of materials to create a new material
with the correct nuclide composition and density. The
:meth:`Material.mix_materials` method takes a list of materials and
a list of their mixing fractions. Mixing fractions can be provided as atomic
fractions, weight fractions, or volume fractions. The fraction type
can be specified by passing 'ao', 'wo', or 'vo' as the third argument, respectively.
For example, assuming the required materials have already been defined, a MOX
material with 3% plutonium oxide by weight could be created using the following:

::

   mox = openmc.Material.mix_materials([uo2, puo2], [0.97, 0.03], 'wo')

It should be noted that, if mixing fractions are specifed as atomic or weight
fractions, the supplied fractions should sum to one. If the fractions are specified
as volume fractions, and the sum of the fractions is less than one, then the remaining
fraction is set as void material.

.. warning:: Materials with :math:`S(\alpha,\beta)` thermal scattering data
             cannot be used in :meth:`Material.mix_materials`. However, thermal
             scattering data can be added to a material created by
             :meth:`Material.mix_materials`.

--------------------
Material Collections
--------------------

The :ref:`scripts_openmc` executable expects to find a ``materials.xml`` file
when it is run. To create this file, one needs to instantiate the
:class:`openmc.Materials` class and add materials to it. The :class:`Materials`
class acts like a list (in fact, it is a subclass of Python's built-in
:class:`list` class), so materials can be added by passing a list to the
constructor, using methods like ``append()``, or through the operator
``+=``. Once materials have been added to the collection, it can be exported
using the :meth:`Materials.export_to_xml` method.

::

   materials = openmc.Materials()
   materials.append(water)
   materials += [uo2, zircaloy]
   materials.export_to_xml()

   # This is equivalent
   materials = openmc.Materials([water, uo2, zircaloy])
   materials.export_to_xml()

Cross Sections
--------------

OpenMC uses a file called :ref:`cross_sections.xml <io_cross_sections>` to
indicate where cross section data can be found on the filesystem. This file
serves the same role that ``xsdir`` does for MCNP_ or ``xsdata`` does for
Serpent. Information on how to generate a cross section listing file can be
found in :ref:`create_xs_library`. Once you have a cross sections file that has
been generated, you can tell OpenMC to use this file either by setting
:attr:`Materials.cross_sections` or by setting the
:envvar:`OPENMC_CROSS_SECTIONS` environment variable to the path of the
``cross_sections.xml`` file. The former approach would look like::

   materials.cross_sections = '/path/to/cross_sections.xml'

.. _MCNP: https://mcnp.lanl.gov/
