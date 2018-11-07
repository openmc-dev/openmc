----------------------------------------------------------
:mod:`openmc.mgxs` -- Multi-Group Cross Section Generation
----------------------------------------------------------

Energy Groups
-------------

Module Variables
++++++++++++++++

.. autodata:: openmc.mgxs.GROUP_STRUCTURES
    :annotation:

Classes
+++++++

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    openmc.mgxs.EnergyGroups

Multi-group Cross Sections
--------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    openmc.mgxs.MGXS
    openmc.mgxs.AbsorptionXS
    openmc.mgxs.CaptureXS
    openmc.mgxs.Chi
    openmc.mgxs.FissionXS
    openmc.mgxs.InverseVelocity
    openmc.mgxs.KappaFissionXS
    openmc.mgxs.MultiplicityMatrixXS
    openmc.mgxs.NuFissionMatrixXS
    openmc.mgxs.ScatterXS
    openmc.mgxs.ScatterMatrixXS
    openmc.mgxs.ScatterProbabilityMatrix
    openmc.mgxs.TotalXS
    openmc.mgxs.TransportXS

Multi-delayed-group Cross Sections
----------------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    openmc.mgxs.MDGXS
    openmc.mgxs.ChiDelayed
    openmc.mgxs.DelayedNuFissionXS
    openmc.mgxs.DelayedNuFissionMatrixXS
    openmc.mgxs.Beta
    openmc.mgxs.DecayRate

Multi-group Cross Section Libraries
-----------------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    openmc.mgxs.Library
