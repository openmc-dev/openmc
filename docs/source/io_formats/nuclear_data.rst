.. _io_nuclear_data:

========================
Nuclear Data File Format
========================

---------------------
Incident Neutron Data
---------------------


**/<nuclide name>/**

:Attributes: - **Z** (*int*) -- Atomic number
             - **A** (*int*) -- Mass number. For a natural element, A=0 is given.
             - **metastable** (*int*) -- Metastable state (0=ground, 1=first
               excited, etc.)
             - **atomic_weight_ratio** (*double*) -- Mass in units of neutron masses
             - **temperature** (*double*) -- Temperature in MeV
             - **n_reaction** (*int*) -- Number of reactions

:Datasets: - **energy** (*double[]*) -- Energy points at which cross sections are tabulated

**/<nuclide name>/reactions/reaction_<mt>/**

:Attributes: - **mt** (*int*) -- ENDF MT reaction number
             - **label** (*char[]*) -- Name of the reaction
             - **Q_value** (*double*) -- Q value in MeV
             - **threshold_idx** (*int*) -- Index on the energy grid that the
               reaction threshold corresponds to
             - **center_of_mass** (*int*) -- Whether the reference frame for
               scattering is center-of-mass (1) or laboratory (0)
             - **n_product** (*int*) -- Number of reaction products

:Datasets: - **xs** (*double[]*) -- Cross section values tabulated against the nuclide energy grid

**/<nuclide name>/reactions/reaction_<mt>/product_<j>/**

   Reaction product data is described in :ref:`product`.

**/<nuclide name>/urr**

:Attributes: - **interpolation** (*int*) -- interpolation scheme
             - **inelastic** (*int*) -- flag indicating inelastic scattering
             - **other_absorb** (*int*) -- flag indicating other absorption
             - **factors** (*int*) -- flag indicating whether tables are
               absolute or multipliers

:Datasets: - **energy** (*double[]*) -- Energy at which probability tables exist
           - **table** (*double[][][]*) -- Probability tables

**/<nuclide name>/total_nu/**

   This special product is used to define the total number of neutrons produced
   from fission. It is formatted as a reaction product, described in
   :ref:`product`.

-------------------------------
Thermal Neutron Scattering Data
-------------------------------

**/<thermal name>/**

:Attributes: - **atomic_weight_ratio** (*double*) -- Mass in units of neutron masses
             - **temperature** (*double*) -- Temperature in MeV
             - **zaids** (*int[]*) -- ZAID identifiers for which the thermal
               scattering data applies to

**/<thermal name>/elastic/**

:Datasets: - **xs** (:ref:`tabulated <1d_tabulated>`) -- Thermal inelastic
             scattering cross section
           - **mu_out** (*double[][]*) -- Distribution of outgoing energies
             and angles for coherent elastic scattering

**/<thermal name>/inelastic/**

:Attributes:
             - **secondary_mode** (*char[]*) -- Indicates how the inelastic
               outgoing angle-energy distributions are represented ('equal',
               'skewed', or 'continuous').

:Datasets: - **xs** (:ref:`tabulated <1d_tabulated>`) -- Thermal inelastic
             scattering cross section
           - **energy_out** (*double[][]*) -- Distribution of outgoing
             energies for each incoming energy. Only present if secondary mode
             is not continuous.
           - **mu_out** (*double[][][]*) -- Distribution of scattering cosines
             for each pair of incoming and outgoing energies. Only present if
             secondary mode is not continuous.

If the secondary mode is continuous, the outgoing energy-angle distribution is
given as a :ref:`correlated angle-energy distribution
<correlated_angle_energy>`.

.. _product:

-----------------
Reaction Products
-----------------

:Object type: Group
:Attributes: - **particle** (*char[]*) -- Type of particle
             - **emission_mode** (*char[]*) -- Emission mode (prompt, delayed,
               total)
             - **decay_rate** (*double*) -- Rate of decay in inverse seconds
             - **n_distribution** (*int*) -- Number of angle/energy
               distributions
:Datasets:
           - **yield** (:ref:`function <1d_functions>`) -- Energy-dependent
             yield of the product.

:Groups:
         - **distribution_<k>** -- Formats for angle-energy distributions are
           detailed in :ref:`angle_energy`. When multiple angle-energy
           distributions occur, one dataset also may appear for each
           distribution:

           :Datasets:
                      - **applicability** (:ref:`function <1d_functions>`) --
                        Probability of selecting this distribution as a function
                        of incident energy

.. _1d_functions:

-------------------------
One-dimensional Functions
-------------------------

Scalar
------

:Object type: Dataset
:Datatype: *double*
:Attributes: - **type** (*char[]*) -- 'constant'

.. _1d_tabulated:

Tabulated
---------

:Object type: Dataset
:Datatype: *double[2][]*
:Description: x-values are listed first followed by corresponding y-values
:Attributes: - **type** (*char[]*) -- 'tabulated'
             - **breakpoints** (*int[]*) -- Region breakpoints
             - **interpolation** (*int[]*) -- Region interpolation codes

Polynomial
----------

:Object type: Dataset
:Datatype: *double[]*
:Description: Polynomial coefficients listed in order of increasing power
:Attributes: - **type** (*char[]*) -- 'polynomial'

Coherent elastic scattering
---------------------------

:Object type: Dataset
:Datatype: *double[2][]*
:Description: The first row lists Bragg edges and the second row lists structure
              factor cumulative sums.
:Attributes: - **type** (*char[]*) -- 'bragg'

.. _angle_energy:

--------------------------
Angle-Energy Distributions
--------------------------

Uncorrelated Angle-Energy
-------------------------

:Object type: Group
:Attributes: - **type** (*char[]*) -- 'uncorrelated'
:Datasets: - **angle/energy** (*double[]*) -- energies at which angle distributions exist
           - **angle/mu** (*double[3][]*) -- tabulated angular distributions for
             each energy. The first row gives :math:`\mu` values, the second row
             gives the probability density, and the third row gives the
             cumulative distribution.

             :Attributes: - **offsets** (*int[]*) -- indices indicating where
                            each angular distribution starts
                          - **interpolation** (*int[]*) -- interpolation code
                            for each angular distribution

:Groups: - **energy/** (:ref:`energy distribution <energy_distribution>`)

.. _correlated_angle_energy:

Correlated Angle-Energy
-----------------------

:Object type: Group
:Attributes: - **type** (*char[]*) -- 'correlated'
:Datasets: - **energy** (*double[]*) -- Incoming energies at which distributions exist

             :Attributes:
                          - **interpolation** (*double[2][]*) -- Breakpoints and
                            interpolation codes for incoming energy regions

           - **energy_out** (*double[5][]*) -- Distribution of outgoing energies
             corresponding to each incoming energy. The distributions are
             flattened into a single array; the start of a given distribution
             can be determined using the ``offsets`` attribute. The first row
             gives outgoing energies, the second row gives the probability
             density, the third row gives the cumulative distribution, the
             fourth row gives interpolation codes for angular distributions, and
             the fifth row gives offsets for angular distributions.

             :Attributes: - **offsets** (*double[]*) -- Offset for each
                            distribution
                          - **interpolation** (*int[]*) -- Interpolation code
                            for each distribution
                          - **n_discrete_lines** (*int[]*) -- Number of discrete
                            lines in each distribution

           - **mu** (*double[3][]*) -- Distribution of angular cosines
             corresponding to each pair of incoming and outgoing energies. The
             distributions are flattened into a single array; the start of a
             given distribution can be determined using offsets in the fifth row
             of the ``energy_out`` dataset. The first row gives angular cosines,
             the second row gives the probability density, and the third row
             gives the cumulative distribution.

Kalbach-Mann
------------

:Object type: Group
:Attributes: - **type** (*char[]*) -- 'kalbach-mann'
:Datasets: - **energy** (*double[]*) -- Incoming energies at which distributions exist

             :Attributes:
                          - **interpolation** (*double[2][]*) -- Breakpoints and
                            interpolation codes for incoming energy regions

           - **distribution** (*double[5][]*) -- Distribution of outgoing
             energies and angles corresponding to each incoming energy. The
             distributions are flattened into a single array; the start of a
             given distribution can be determined using the ``offsets``
             attribute. The first row gives outgoing energies, the second row
             gives the probability density, the third row gives the cumulative
             distribution, the fourth row gives Kalbach-Mann precompound
             factors, and the fifth row gives Kalbach-Mann angular distribution
             slopes.

             :Attributes: - **offsets** (*double[]*) -- Offset for each
                            distribution
                          - **interpolation** (*int[]*) -- Interpolation code
                            for each distribution
                          - **n_discrete_lines** (*int[]*) -- Number of discrete
                            lines in each distribution

N-Body Phase Space
------------------

:Object type: Group
:Attributes: - **type** (*char[]*) -- 'nbody'
             - **total_mass** (*double*) -- Total mass of product particles
             - **n_particles** (*int*) -- Number of product particles
             - **atomic_weight_ratio** (*double*) -- Atomic weight ratio of the
               target nuclide in neutron masses
             - **q_value** (*double*) -- Q value for the reaction in MeV

.. _energy_distribution:

--------------------
Energy Distributions
--------------------

Maxwell
-------

:Object type: Group
:Attributes: - **type** (*char[]*) -- 'maxwell'
             - **u** (*double*) -- Restriction energy in MeV
:Datasets:
           - **theta** (:ref:`tabulated <1d_tabulated>`) -- Maxwellian
             temperature as a function of energy

Evaporation
-----------

:Object type: Group
:Attributes: - **type** (*char[]*) -- 'evaporation'
             - **u** (*double*) -- Restriction energy in MeV
:Datasets:
           - **theta** (:ref:`tabulated <1d_tabulated>`) -- Evaporation
             temperature as a function of energy

Watt Fission Spectrum
---------------------

:Object type: Group
:Attributes: - **type** (*char[]*) -- 'watt'
             - **u** (*double*) -- Restriction energy in MeV
:Datasets: - **a** (:ref:`tabulated <1d_tabulated>`) -- Watt parameter :math:`a`
             as a function of incident energy
           - **b** (:ref:`tabulated <1d_tabulated>`) -- Watt parameter :math:`b`
             as a function of incident energy

Madland-Nix
-----------

:Object type: Group
:Attributes: - **type** (*char[]*) -- 'watt'
             - **efl** (*double*) -- Average energy of light fragment in eV
             - **efh** (*double*) -- Average energy of heavy fragment in eV

Discrete Photon
---------------

:Object type: Group
:Attributes: - **type** (*char[]*) -- 'discrete_photon'
             - **primary_flag** (*int*) -- Whether photon is a primary
             - **energy** (*double*) -- Photon energy in MeV
             - **atomic_weight_ratio** (*double*) -- Atomic weight ratio of
               target nuclide in neutron masses

Level Inelastic
---------------

:Object type: Group
:Attributes: - **type** (*char[]*) -- 'level'
             - **threshold** (*double*) -- Energy threshold in the laboratory
               system in MeV
             - **mass_ratio** (*double*) -- :math:`(A/(A + 1))^2`

Continuous Tabular
------------------

:Object type: Group
:Attributes: - **type** (*char[]*) -- 'continuous'
:Datasets: - **energy** (*double[]*) -- Incoming energies at which distributions exist

             :Attributes:
                          - **interpolation** (*double[2][]*) -- Breakpoints and
                            interpolation codes for incoming energy regions

           - **distribution** (*double[3][]*) -- Distribution of outgoing
             energies corresponding to each incoming energy. The distributions
             are flattened into a single array; the start of a given
             distribution can be determined using the ``offsets`` attribute. The
             first row gives outgoing energies, the second row gives the
             probability density, and the third row gives the cumulative
             distribution.

             :Attributes: - **offsets** (*double[]*) -- Offset for each
                            distribution
                          - **interpolation** (*int[]*) -- Interpolation code
                            for each distribution
                          - **n_discrete_lines** (*int[]*) -- Number of discrete
                            lines in each distribution
