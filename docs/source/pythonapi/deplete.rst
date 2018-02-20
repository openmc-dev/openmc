.. _pythonapi_deplete:

----------------------------------
:mod:`openmc.deplete` -- Depletion
----------------------------------

Integrators
-----------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   openmc.deplete.integrator.predictor
   openmc.deplete.integrator.cecm

Integrator Helper Functions
---------------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   openmc.deplete.integrator.CRAM16
   openmc.deplete.integrator.CRAM48
   openmc.deplete.integrator.save_results

Metaclasses
-----------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.deplete.TransportOperator

OpenMC Classes
--------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.deplete.Operator
   openmc.deplete.OperatorResult

Data Classes
------------
.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.deplete.AtomNumber
   openmc.deplete.Chain
   openmc.deplete.Nuclide
   openmc.deplete.ReactionRates
   openmc.deplete.Results
