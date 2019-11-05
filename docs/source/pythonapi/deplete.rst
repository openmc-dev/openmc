.. _pythonapi_deplete:

.. module:: openmc.deplete

----------------------------------
:mod:`openmc.deplete` -- Depletion
----------------------------------

Primary API
-----------

The two primary requirements to perform depletion with :mod:`openmc.deplete`
are:

    1) A transport operator
    2) A time-integration scheme

The former is responsible for executing a transport code, like OpenMC,
and retaining important information required for depletion. The most common examples
are reaction rates and power normalization data. The latter is responsible for
projecting reaction rates and compositions forward in calendar time across
some step size :math:`\Delta t`, and obtaining new compositions given a power
or power density. The :class:`Operator` is provided to handle communicating with
OpenMC. Several classes are provided that implement different time-integration
algorithms for depletion calculations, which are described in detail in Colin
Josey's thesis, `Development and analysis of high order neutron
transport-depletion coupling algorithms <http://hdl.handle.net/1721.1/113721>`_.

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template:  myintegrator.rst

    PredictorIntegrator
    CECMIntegrator
    CELIIntegrator
    CF4Integrator
    EPCRK4Integrator
    LEQIIntegrator
    SICELIIntegrator
    SILEQIIntegrator

Each of these classes expects a "transport operator" to be passed. An operator
specific to OpenMC is available using the following class:

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: mycallable.rst

   Operator

The :class:`Operator` must also have some knowledge of how nuclides transmute
and decay. This is handled by the :class:`Chain`.

Minimal Example
---------------

A minimal example for performing depletion would be:

.. code::

    >>> import openmc
    >>> import openmc.deplete
    >>> geometry = openmc.Geometry.from_xml()
    >>> settings = openmc.Settings.from_xml()

    # Representation of a depletion chain
    >>> chain_file = "chain_casl.xml"
    >>> operator = openmc.deplete.Operator(
    ...     geometry, settings, chain_file)

    # Set up 5 time steps of one day each
    >>> dt = [24 * 60 * 60] * 5
    >>> power = 1e6  # constant power of 1 MW

    # Deplete using mid-point predictor-corrector
    >>> cecm = openmc.deplete.CECMIntegrator(
    ...     operator, dt, power)
    >>> cecm.integrate()

Internal Classes and Functions
------------------------------

When running in parallel using `mpi4py
<https://mpi4py.readthedocs.io/en/stable/>`_, the MPI intercommunicator used can
be changed by modifying the following module variable. If it is not explicitly
modified, it defaults to ``mpi4py.MPI.COMM_WORLD``.

.. data:: comm

   MPI intercommunicator used to call OpenMC library

   :type: mpi4py.MPI.Comm

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
   FissionYieldDistribution
   FissionYield

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

The following class and functions are used to solve the depletion equations,
with :func:`cram.CRAM48` being the default.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myintegrator.rst

   cram.IPFCramSolver

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   cram.CRAM16
   cram.CRAM48
   cram.deplete
   cram.timed_deplete

The following classes are used to help the :class:`openmc.deplete.Operator`
compute quantities like effective fission yields, reaction rates, and
total system energy.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   helpers.AveragedFissionYieldHelper
   helpers.ChainFissionHelper
   helpers.ConstantFissionYieldHelper
   helpers.DirectReactionRateHelper
   helpers.EnergyScoreHelper
   helpers.FissionYieldCutoffHelper


Abstract Base Classes
---------------------

A good starting point for extending capabilities in :mod:`openmc.deplete` is
to examine the following abstract base classes. Custom classes can
inherit from :class:`abc.TransportOperator` to implement alternative
schemes for collecting reaction rates and other data from a transport code
prior to depleting materials

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: mycallable.rst

   abc.TransportOperator

The following classes are abstract classes used to pass information from
OpenMC simulations back on to the :class:`abc.TransportOperator`

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   abc.EnergyHelper
   abc.FissionYieldHelper
   abc.ReactionRateHelper
   abc.TalliedFissionYieldHelper

Custom integrators or depletion solvers can be developed by subclassing from
the following abstract base classes:

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myintegrator.rst

   abc.Integrator
   abc.SIIntegrator
   abc.DepSystemSolver
