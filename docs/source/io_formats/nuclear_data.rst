.. _io_nuclear_data:

=========================
Nuclear Data File Formats
=========================

---------------------
Incident Neutron Data
---------------------

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file
             - **version** (*int[2]*) -- Major and minor version of the data

**/<nuclide name>/**

:Attributes: - **Z** (*int*) -- Atomic number
             - **A** (*int*) -- Mass number. For a natural element, A=0 is given.
             - **metastable** (*int*) -- Metastable state (0=ground, 1=first
               excited, etc.)
             - **atomic_weight_ratio** (*double*) -- Mass in units of neutron masses
             - **n_reaction** (*int*) -- Number of reactions

:Datasets:
           - **energy** (*double[]*) -- Energies in [eV] at which cross sections
             are tabulated

**/<nuclide name>/kTs/**

<TTT>K is the temperature in Kelvin, rounded to the nearest integer, of the
temperature-dependent data set.  For example, the data set corresponding to
300 Kelvin would be located at `300K`.

:Datasets:
           - **<TTT>K** (*double*) -- kT values in [eV] for each temperature
             TTT (in Kelvin)

**/<nuclide name>/reactions/reaction_<mt>/**

:Attributes: - **mt** (*int*) -- ENDF MT reaction number
             - **label** (*char[]*) -- Name of the reaction
             - **Q_value** (*double*) -- Q value in eV
             - **center_of_mass** (*int*) -- Whether the reference frame for
               scattering is center-of-mass (1) or laboratory (0)
             - **n_product** (*int*) -- Number of reaction products
             - **redundant** (*int*) -- Whether reaction is redundant

**/<nuclide name>/reactions/reaction_<mt>/<TTT>K/**

<TTT>K is the temperature in Kelvin, rounded to the nearest integer, of the
temperature-dependent data set.  For example, the data set corresponding to
300 Kelvin would be located at `300K`.

:Datasets:
           - **xs** (*double[]*) -- Cross section values tabulated against the
             nuclide energy grid for temperature TTT (in Kelvin)

             :Attributes:
                          - **threshold_idx** (*int*) -- Index on the energy
                            grid that the reaction threshold corresponds to for
                            temperature TTT (in Kelvin)

**/<nuclide name>/reactions/reaction_<mt>/product_<j>/**

   Reaction product data is described in :ref:`product`.

**/<nuclide name>/urr/<TTT>K/**

<TTT>K is the temperature in Kelvin, rounded to the nearest integer, of the
temperature-dependent data set.  For example, the data set corresponding to
300 Kelvin would be located at `300K`.

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

**/<nuclide name>/fission_energy_release/**

:Datasets: - **fragments** (:ref:`function <1d_functions>`) -- Energy
             released in the form of fragments as a function of incident
             neutron energy.
           - **prompt_neutrons** (:ref:`function <1d_functions>`) -- Energy
             released in the form of prompt neutrons as a function of incident
             neutron energy.
           - **delayed_neutrons** (:ref:`function <1d_functions>`) -- Energy
             released in the form of delayed neutrons as a function of incident
             neutron energy.
           - **prompt_photons** (:ref:`function <1d_functions>`) -- Energy
             released in the form of prompt photons as a function of incident
             neutron energy.
           - **delayed_photons** (:ref:`function <1d_functions>`) -- Energy
             released in the form of delayed photons as a function of incident
             neutron energy.
           - **betas** (:ref:`function <1d_functions>`) -- Energy released in
             the form of betas as a function of incident neutron energy.
           - **neutrinos** (:ref:`function <1d_functions>`) -- Energy released
             in the form of neutrinos as a function of incident neutron energy.
           - **q_prompt** (:ref:`function <1d_functions>`) -- The prompt fission
             Q-value (fragments + prompt neutrons + prompt photons - incident
             energy)
           - **q_recoverable** (:ref:`function <1d_functions>`) -- The
             recoverable fission Q-value (Q_prompt + delayed neutrons + delayed
             photons + betas)

--------------------
Incident Photon Data
--------------------

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file
             - **version** (*int[2]*) -- Major and minor version of the data

**/<element>/**

:Attributes: - **Z** (*int*) -- Atomic number

:Datasets:
           - **energy** (*double[]*) -- Energies in [eV] at which cross sections
             are tabulated

**/<element>/bremsstrahlung/**

:Attributes: - **I** (*double*) -- Mean excitation energy in [eV]

:Datasets: - **electron_energy** (*double[]*) -- Incident electron energy in [eV]
           - **photon_energy** (*double[]*) -- Outgoing photon energy as
             fraction of incident electron energy
           - **dcs** (*double[][]*) -- Bremsstrahlung differential cross section
             at each incident energy in [mb/eV]
           - **ionization_energy** (*double[]*) -- Ionization potential of each
             subshell in [eV]
           - **num_electrons** (*int[]*) -- Number of electrons per subshell,
             with conduction electrons indicated by a negative value

**/<element>/coherent/**

:Datasets: - **xs** (*double[]*) -- Coherent scattering cross section in [b]
           - **integrated_scattering_factor** (:ref:`tabulated <1d_tabulated>`)
             -- Integrated coherent scattering form factor
           - **anomalous_real** (:ref:`tabulated <1d_tabulated>`) -- Real part
             of the anomalous scattering factor
           - **anomalous_imag** (:ref:`tabulated <1d_tabulated>`) -- Imaginary
             part of the anomalous scattering factor

**/<element>/compton_profiles/**

:Datasets: - **binding_energy** (*double[]*) -- Binding energy for each subshell in [eV]
           - **num_electrons** (*double[]*) -- Number of electrons in each subshell
           - **pz** (*double[]*) -- Projection of the electron momentum on the
             scattering vector in units of :math:`me^2 / \hbar` where :math:`m`
             is the electron rest mass and :math:`e` is the electron charge
           - **J** (*double[][]*) -- Compton profile for each subshell in units
             of :math:`\hbar / (me^2)`

**/<element>/heating/**

:Datasets: - **xs** (*double[]*) -- Total heating cross section in [b-eV]

**/<element>/incoherent/**

:Datasets: - **xs** (*double[]*) -- Incoherent scattering cross section in [b]
           - **scattering_factor** (:ref:`tabulated <1d_tabulated>`) --

**/<element>/pair_production_electron/**

:Datasets: - **xs** (*double[]*) -- Pair production (electron field) cross section in [b]

**/<element>/pair_production_nuclear/**

:Datasets: - **xs** (*double[]*) -- Pair production (nuclear field) cross section in [b]

**/<element>/photoelectric/**

:Datasets: - **xs** (*double[]*) -- Total photoionization cross section in [b]

**/<element>/subshells/**

:Attributes: - **designators** (*char[][]*) -- Designator for each shell, e.g. 'M2'

**/<element>/subshells/<designator>/**

:Attributes: - **binding_energy** (*double*) -- Binding energy of the subshell in [eV]
             - **num_electrons** (*double*) -- Number of electrons in the subshell

:Datasets: - **transitions** (*double[][]*) -- Atomic relaxation data
           - **xs** (*double[]*) -- Photoionization cross section for subshell
             in [b] tabulated against the main energy grid

             :Attributes:
                          - **threshold_idx** (*int*) -- Index on the energy
                            grid of the reaction threshold

-------------------------------
Thermal Neutron Scattering Data
-------------------------------

**/**

:Attributes:
             - **version** (*int[2]*) -- Major and minor version of the data

**/<thermal name>/**

:Attributes: - **atomic_weight_ratio** (*double*) -- Mass in units of neutron masses
             - **energy_max** (*double*) -- Maximum energy in [eV]
             - **nuclides** (*char[][]*) -- Names of nuclides for which the
               thermal scattering data applies to

**/<thermal name>/kTs/**

<TTT>K is the temperature in Kelvin, rounded to the nearest integer, of the
temperature-dependent data set.  For example, the data set corresponding to
300 Kelvin would be located at `300K`.

:Datasets:
           - **<TTT>K** (*double*) -- kT values (in eV) for each temperature
             TTT (in Kelvin)

**/<thermal name>/elastic/<TTT>K/**

<TTT>K is the temperature in Kelvin, rounded to the nearest integer, of the
temperature-dependent data set.  For example, the data set corresponding to
300 Kelvin would be located at `300K`.

:Datasets:
           - **xs** (:ref:`function <1d_functions>`) -- Thermal elastic
             scattering cross section for temperature TTT (in Kelvin)

:Groups:
         - **distribution** -- Format for angle-energy distributions are
           detailed in :ref:`angle_energy`.

**/<thermal name>/inelastic/<TTT>K/**

<TTT>K is the temperature in Kelvin, rounded to the nearest integer, of the
temperature-dependent data set.  For example, the data set corresponding to
300 Kelvin would be located at `300K`.

:Datasets:
           - **xs** (:ref:`function <1d_functions>`) -- Thermal inelastic
             scattering cross section for temperature TTT (in Kelvin)

:Groups:
         - **distribution** -- Format for angle-energy distributions are
           detailed in :ref:`angle_energy`.

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
:Attributes: - **type** (*char[]*) -- 'Tabulated1D'
             - **breakpoints** (*int[]*) -- Region breakpoints
             - **interpolation** (*int[]*) -- Region interpolation codes

.. _1d_polynomial:

Polynomial
----------

:Object type: Dataset
:Datatype: *double[]*
:Description: Polynomial coefficients listed in order of increasing power
:Attributes: - **type** (*char[]*) -- 'Polynomial'

Coherent elastic scattering
---------------------------

:Object type: Dataset
:Datatype: *double[2][]*
:Description: The first row lists Bragg edges and the second row lists structure
              factor cumulative sums.
:Attributes: - **type** (*char[]*) -- 'CoherentElastic'

Incoherent elastic scattering
-----------------------------

:Object type: Dataset
:Datatype: *double[2]*
:Description: The first value is the characteristic bound cross section in [b]
              and the second value is the Debye-Waller integral in
              [eV\ :math:`^{-1}`].
:Attributes: - **type** (*char[]*) -- 'IncoherentElastic'

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
             - **q_value** (*double*) -- Q value for the reaction in eV

Coherent Elastic
----------------

This angle-energy distribution is used specifically for coherent elastic thermal
neutron scattering.

:Object type: Group
:Attributes: - **type** (*char[]*) -- "coherent_elastic"
:Hard link: - **xs** -- Link to the coherent elastic scattering cross section

Incoherent Elastic
------------------

This angle-energy distribution is used specifically for incoherent elastic
thermal neutron scattering (derived from an ENDF file directly).

:Object type: Group
:Attributes: - **type** (*char[]*) -- "incoherent_elastic"
:Datasets:
           - **debye_waller** (*double*) -- Debye-Waller integral in
             [eV\ :math:`^{-1}`]

Incoherent Elastic (Discrete)
-----------------------------

This angle-energy distribution is used for discretized incoherent elastic
thermal neutron scattering distributions that are present in ACE files.

:Object type: Group
:Attributes: - **type** (*char[]*) -- "incoherent_elastic_discrete"
:Datasets:
           - **mu_out** (*double[][]*) -- Equiprobable discrete outgoing
             angles for each incident neutron energy tabulated

Incoherent Inelastic
--------------------

This angle-energy distribution is used specifically for (continuous) incoherent
inelastic thermal neutron scattering.

:Object type: Group
:Attributes: - **type** (*char[]*) -- "incoherent_inelastic"
:Datasets: The datasets for this angle-energy distribution are the same as for
           :ref:`correlated angle-energy distributions
           <correlated_angle_energy>`.

Incoherent Inelastic (Discrete)
-------------------------------

This angle-energy distribution is used specifically for incoherent inelastic
thermal neutron scattering where the distributions have been discretized into
equiprobable bins.

:Object type: Group
:Attributes: - **type** (*char[]*) -- "incoherent_inelastic_discrete"
:Datasets: - **energy_out** (*double[][]*) -- Distribution of outgoing
             energies for each incoming energy.
           - **mu_out** (*double[][][]*) -- Distribution of scattering cosines
             for each pair of incoming and outgoing energies.
           - **skewed** (*int8_t*) -- Whether discrete angles are equi-probable
             (0) or have a skewed distribution (1).

.. _energy_distribution:

--------------------
Energy Distributions
--------------------

Maxwell
-------

:Object type: Group
:Attributes: - **type** (*char[]*) -- 'maxwell'
             - **u** (*double*) -- Restriction energy in eV
:Datasets:
           - **theta** (:ref:`tabulated <1d_tabulated>`) -- Maxwellian
             temperature as a function of energy

Evaporation
-----------

:Object type: Group
:Attributes: - **type** (*char[]*) -- 'evaporation'
             - **u** (*double*) -- Restriction energy in eV
:Datasets:
           - **theta** (:ref:`tabulated <1d_tabulated>`) -- Evaporation
             temperature as a function of energy

Watt Fission Spectrum
---------------------

:Object type: Group
:Attributes: - **type** (*char[]*) -- 'watt'
             - **u** (*double*) -- Restriction energy in eV
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
             - **energy** (*double*) -- Photon energy in eV
             - **atomic_weight_ratio** (*double*) -- Atomic weight ratio of
               target nuclide in neutron masses

Level Inelastic
---------------

:Object type: Group
:Attributes: - **type** (*char[]*) -- 'level'
             - **threshold** (*double*) -- Energy threshold in the laboratory
               system in eV
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
