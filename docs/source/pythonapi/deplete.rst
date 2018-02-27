.. _pythonapi_deplete:

----------------------------------
:mod:`openmc.deplete` -- Depletion
----------------------------------

.. module:: openmc.deplete

Two functions are provided that implement different time-integration algorithms
for depletion calculations.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   integrator.predictor
   integrator.cecm

Each of these functions expects a "transport operator" to be passed. An operator
specific to OpenMC is available using the following class:

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   Operator

When running in parallel using `mpi4py <http://mpi4py.scipy.org>`_, the MPI
intercommunicator used can be changed by modifying the following module
variable. If it is not explicitly modified, it defaults to
``mpi4py.MPI.COMM_WORLD``.

.. data:: comm

   MPI intercommunicator used to call OpenMC library

   :type: mpi4py.MPI.Comm

Internal Classes and Functions
------------------------------

During a depletion calculation, the depletion chain, reaction rates, and number
densities are managed through a series of internal classes that are not normally
visible to a user. However, should you find yourself wondering about these
classes (e.g., if you want to know what decay modes or reactions are present in
a depletion chain), they are documented here. The following classes store data
for a depletion chain:

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   Chain
   DecayTuple
   Nuclide
   ReactionTuple

The following classes are used during a depletion simulation and store auxiliary
data, such as number densities and reaction rates for each material.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   AtomNumber
   OperatorResult
   ReactionRates
   Results
   ResultsList
   TransportOperator

Each of the integrator functions also relies on a number of "helper" functions
as follows:

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   integrator.CRAM16
   integrator.CRAM48
