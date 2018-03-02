--------------------------------------------
:mod:`openmc.data` -- Nuclear Data Interface
--------------------------------------------

Core Classes
------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    openmc.data.IncidentNeutron
    openmc.data.Reaction
    openmc.data.Product
    openmc.data.Tabulated1D
    openmc.data.FissionEnergyRelease
    openmc.data.ThermalScattering
    openmc.data.CoherentElastic
    openmc.data.FissionEnergyRelease
    openmc.data.DataLibrary
    openmc.data.Decay
    openmc.data.FissionProductYields
    openmc.data.WindowedMultipole

Core Functions
--------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myfunction.rst

    openmc.data.atomic_mass
    openmc.data.gnd_name
    openmc.data.linearize
    openmc.data.thin
    openmc.data.water_density
    openmc.data.write_compact_458_library
    openmc.data.zam

Angle-Energy Distributions
--------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    openmc.data.AngleEnergy
    openmc.data.KalbachMann
    openmc.data.CorrelatedAngleEnergy
    openmc.data.UncorrelatedAngleEnergy
    openmc.data.NBodyPhaseSpace
    openmc.data.LaboratoryAngleEnergy
    openmc.data.AngleDistribution
    openmc.data.EnergyDistribution
    openmc.data.ArbitraryTabulated
    openmc.data.GeneralEvaporation
    openmc.data.MaxwellEnergy
    openmc.data.Evaporation
    openmc.data.WattEnergy
    openmc.data.MadlandNix
    openmc.data.DiscretePhoton
    openmc.data.LevelInelastic
    openmc.data.ContinuousTabular

Resonance Data
--------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    openmc.data.Resonances
    openmc.data.ResonanceRange
    openmc.data.SingleLevelBreitWigner
    openmc.data.MultiLevelBreitWigner
    openmc.data.ReichMoore
    openmc.data.RMatrixLimited
    openmc.data.ParticlePair
    openmc.data.SpinGroup
    openmc.data.Unresolved

ACE Format
----------

Classes
+++++++

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    openmc.data.ace.Library
    openmc.data.ace.Table

Functions
+++++++++

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myfunction.rst

    openmc.data.ace.ascii_to_binary

ENDF Format
-----------

Classes
+++++++

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    openmc.data.endf.Evaluation

Functions
+++++++++

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myfunction.rst

    openmc.data.endf.float_endf
    openmc.data.endf.get_cont_record
    openmc.data.endf.get_evaluations
    openmc.data.endf.get_head_record
    openmc.data.endf.get_tab1_record
    openmc.data.endf.get_tab2_record
    openmc.data.endf.get_text_record

NJOY Interface
--------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myfunction.rst

    openmc.data.njoy.run
    openmc.data.njoy.make_pendf
    openmc.data.njoy.make_ace
    openmc.data.njoy.make_ace_thermal
