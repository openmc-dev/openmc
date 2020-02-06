--------------------------------------------
:mod:`openmc.data` -- Nuclear Data Interface
--------------------------------------------

.. module:: openmc.data

Core Classes
------------

The following classes are used for incident neutron data, decay data, fission
and product yields.

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    IncidentNeutron
    Reaction
    Product
    FissionEnergyRelease
    DataLibrary
    Decay
    FissionProductYields
    WindowedMultipole
    ProbabilityTables

The following classes are used for storing atomic data (incident photon cross
sections, atomic relaxation):

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    IncidentPhoton
    PhotonReaction
    AtomicRelaxation


The following classes are used for storing thermal neutron scattering data:

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    ThermalScattering
    ThermalScatteringReaction
    CoherentElastic
    IncoherentElastic


Core Functions
--------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myfunction.rst

    atomic_mass
    gnd_name
    linearize
    thin
    water_density
    zam

One-dimensional Functions
-------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    Function1D
    Tabulated1D
    Polynomial
    Combination
    Sum
    Regions1D
    ResonancesWithBackground

Angle-Energy Distributions
--------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    AngleEnergy
    KalbachMann
    CorrelatedAngleEnergy
    UncorrelatedAngleEnergy
    NBodyPhaseSpace
    LaboratoryAngleEnergy
    AngleDistribution
    EnergyDistribution
    ArbitraryTabulated
    GeneralEvaporation
    MaxwellEnergy
    Evaporation
    WattEnergy
    MadlandNix
    DiscretePhoton
    LevelInelastic
    ContinuousTabular
    CoherentElasticAE
    IncoherentElasticAE
    IncoherentElasticAEDiscrete
    IncoherentInelasticAEDiscrete

Resonance Data
--------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    Resonances
    ResonanceRange
    SingleLevelBreitWigner
    MultiLevelBreitWigner
    ReichMoore
    RMatrixLimited
    ResonanceCovariances
    ResonanceCovarianceRange
    SingleLevelBreitWignerCovariance
    MultiLevelBreitWignerCovariance
    ReichMooreCovariance
    ParticlePair
    SpinGroup
    Unresolved

ACE Format
----------

Classes
+++++++

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    ace.Library
    ace.Table
    ace.TableType

Functions
+++++++++

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myfunction.rst

    ace.ascii_to_binary
    ace.get_libraries_from_xsdir
    ace.get_libraries_from_xsdata

ENDF Format
-----------

Classes
+++++++

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    endf.Evaluation

Functions
+++++++++

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myfunction.rst

    endf.float_endf
    endf.get_cont_record
    endf.get_evaluations
    endf.get_head_record
    endf.get_tab1_record
    endf.get_tab2_record
    endf.get_text_record

NJOY Interface
--------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myfunction.rst

    njoy.run
    njoy.make_pendf
    njoy.make_ace
    njoy.make_ace_thermal
