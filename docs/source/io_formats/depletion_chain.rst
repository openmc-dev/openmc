.. _io_depletion_chain:

============================
Depletion Chain -- chain.xml
============================

A depletion chain file has a ``<depletion_chain>`` root element with one or more
``<nuclide>`` child elements. The decay, reaction, and fission product data for
each nuclide appears as child elements of ``<nuclide>``.

.. _io_chain_source_metadata:

-----------------------------
``<source_metadata>`` Element
-----------------------------

An optional ``<source_metadata>`` child of ``<depletion_chain>`` records the
ENDF evaluations supplied when constructing the chain. It contains one
``<source>`` element for each distinct library, version, and release in each
component. These source records are separate from the particle-emission
``<source>`` elements inside a ``<nuclide>``.

Each metadata ``<source>`` has the following required attributes:

  :component:
    ``neutron``, ``decay``, or ``fission_yield``.

  :library:
    ENDF library name, such as ``ENDF/B``.

  :version:
    Nonnegative integer version from the ENDF evaluation header.

  :release:
    Nonnegative integer release from the ENDF evaluation header.

For example, a chain constructed from three ENDF/B-VII.1 sublibraries can
contain:

.. code-block:: xml

    <source_metadata>
      <source component="neutron" library="ENDF/B" version="7" release="1"/>
      <source component="decay" library="ENDF/B" version="7" release="1"/>
      <source component="fission_yield" library="ENDF/B" version="7" release="1"/>
    </source_metadata>

A component may have several records when its inputs come from different
libraries or releases. Missing metadata means the source was not recorded;
it does not imply a particular library. Existing chain files without this
element remain supported.

This is construction provenance, not a selection of transport cross sections
or a record of subsequent processing and manual edits. Chain reduction
preserves these original construction records, including sources for nuclides
that may have been removed.

---------------------
``<nuclide>`` Element
---------------------

The ``<nuclide>`` element contains information on the decay modes, reactions,
and fission product yields for a given nuclide in the depletion chain. This
element may have the following attributes:

  :name:
    Name of the nuclide

  :half_life:
    Half-life of the nuclide in [s]

  :decay_modes:
    Number of decay modes present

  :decay_energy:
    Decay energy released in [eV]

  :reactions:
    Number of reactions present

For each decay mode, a :ref:`io_chain_decay` appears as a child of
``<nuclide>``. For each reaction present, a :ref:`io_chain_reaction` appears as
a child of ``<nuclide>``. If the nuclide is fissionable, a :ref:`io_chain_nfy`
appears as well.

.. _io_chain_decay:

-------------------
``<decay>`` Element
-------------------

The ``<decay>`` element represents a single decay mode and has the following
attributes:

  :type:
    The type of the decay, e.g. 'ec/beta+'

  :target:
    The daughter nuclide produced from the decay

  :branching_ratio:
    The branching ratio for this decay mode

.. _io_chain_reaction:

--------------------
``<source>`` Element
--------------------

The ``<source>`` element represents photon and electron sources associated with
the decay of a nuclide and contains information to construct an
:class:`openmc.stats.Univariate` object that represents this emission as an
energy distribution. This element has the following attributes:

  :type:
    The type of :class:`openmc.stats.Univariate` source term.

  :particle:
    The type of particle emitted, e.g., 'photon' or 'electron'

  :parameters:
    The parameters of the source term, e.g., for a
    :class:`openmc.stats.Discrete` source, the energies (in [eV]) at which the
    particles are emitted and their relative intensities in [Bq/atom] (in other
    words, decay constants).

----------------------
``<reaction>`` Element
----------------------

The ``<reaction>`` element represents a single transmutation reaction. This
element has the following attributes:

  :type:
    The type of the reaction, e.g., '(n,gamma)'

  :Q:
    The Q value of the reaction in [eV]

  :target:
    The nuclide produced in the reaction (absent if the type is 'fission')

  :branching_ratio:
    The branching ratio for the reaction

.. _io_chain_nfy:

------------------------------------
``<neutron_fission_yields>`` Element
------------------------------------

The ``<neutron_fission_yields>`` element provides yields of fission products for
fissionable nuclides. Normally, it has the follow sub-elements:

  :energies:
    Energies in [eV] at which yields for products are tabulated

  :fission_yields:

    Fission product yields for a single energy point. This element itself has a
    number of attributes/sub-elements:

      :energy:
        Energy in [eV] at which yields are tabulated

      :products:
        Names of fission products

      :data:
        Independent yields for each fission product

In the event that a nuclide doesn't have any known fission product yields, it is
possible to have that nuclide borrow yields from another nuclide by indicating
the other nuclide in a single `parent` attribute. For example:

.. code-block:: xml

    <neutron_fission_yields parent="U235"/>
