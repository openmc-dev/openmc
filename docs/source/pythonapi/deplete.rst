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

The former is responsible for calculating and retaining important information
required for depletion. The most common examples are reaction rates and power
normalization data. The latter is responsible for projecting reaction rates and
compositions forward in calendar time across some step size :math:`\Delta t`,
and obtaining new compositions given a power or power density. The
:class:`CoupledOperator` class is provided to obtain reaction rates via tallies
through OpenMC's transport solver, and the :class:`IndependentOperator` class is
provided to obtain reaction rates from cross-section data. Several classes are
provided that implement different time-integration algorithms for depletion
calculations, which are described in detail in Colin Josey's thesis,
`Development and analysis of high order neutron transport-depletion coupling
algorithms <https://dspace.mit.edu/handle/1721.1/113721>`_.

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

Each of these classes expects a "transport operator" to be passed. OpenMC
provides the following transport operator classes:

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: mycallable.rst

   CoupledOperator
   IndependentOperator

The :class:`CoupledOperator` and :class:`IndependentOperator` classes must also
have some knowledge of how nuclides transmute and decay. This is handled by the
:class:`Chain` class.

The :class:`IndependentOperator` class requires a set of fluxes and microscopic
cross sections. The following function can be used to generate this information:

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   get_microxs_and_flux

Minimal Example
---------------

A minimal example for performing depletion would be:

.. code::

    >>> import openmc
    >>> import openmc.deplete
    >>> geometry = openmc.Geometry.from_xml()
    >>> settings = openmc.Settings.from_xml()
    >>> model = openmc.model.Model(geometry, settings)

    # Representation of a depletion chain
    >>> chain_file = "chain_casl.xml"
    >>> operator = openmc.deplete.CoupledOperator(
    ...     model, chain_file)

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

The :class:`Chain` class uses information from the following module variable:

.. data:: chain.REACTIONS

   Dictionary that maps transmutation reaction names to information needed when
   a chain is being generated: MT values, the change in atomic/mass numbers
   resulting from the reaction, and what secondaries are produced.

   :type: dict

The following classes are used during a depletion simulation and store auxiliary
data, such as number densities and reaction rates for each material.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   AtomNumber
   MicroXS
   OperatorResult
   ReactionRates
   Results
   StepResult

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
   pool.deplete

.. data:: pool.USE_MULTIPROCESSING

   Boolean switch to enable or disable the use of :mod:`multiprocessing`
   when solving the Bateman equations. The default is to use
   :mod:`multiprocessing`, but can cause the simulation to hang in
   some computing environments, namely due to MPI and networking
   restrictions. Disabling this option will result in only a single
   CPU core being used for depletion.

   :type: bool

.. data:: pool.NUM_PROCESSES

   Number of worker processes used for depletion calculations, which rely on the
   :class:`multiprocessing.pool.Pool` class. If set to ``None`` (default), the
   number returned by :func:`os.cpu_count` is used.

The following classes are used to help the :class:`openmc.deplete.CoupledOperator`
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
   helpers.FluxCollapseHelper

The :class:`openmc.deplete.IndependentOperator` uses inner classes subclassed
from those listed above to perform similar calculations.

The following classes are used to define transfer rates to model continuous
removal or feed of nuclides during depletion.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   transfer_rates.TransferRates

Intermediate Classes
--------------------

Specific implementations of abstract base classes may utilize some of
the same methods and data structures. These methods and data are stored
in intermediate classes.

Methods common to tally-based implementation of :class:`FissionYieldHelper`
are stored in :class:`helpers.TalliedFissionYieldHelper`

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   helpers.TalliedFissionYieldHelper

Methods common to OpenMC-specific implementations of :class:`TransportOperator`
are stored in :class:`openmc_operator.OpenMCOperator`

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: mycallable.rst

   openmc_operator.OpenMCOperator


Abstract Base Classes
---------------------

A good starting point for extending capabilities in :mod:`openmc.deplete` is
to examine the following abstract base classes. Custom classes can
inherit from :class:`abc.TransportOperator` to implement alternative
schemes for collecting reaction rates and other data prior to depleting
materials

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: mycallable.rst

   abc.TransportOperator

The following classes are abstract classes used to pass information from
transport simulations (in the case of transport-coupled depletion) or to
simply calculate these quantities directly (in the case of
transport-independent depletion) back on to the :class:`abc.TransportOperator`

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   abc.NormalizationHelper
   abc.FissionYieldHelper
   abc.ReactionRateHelper

Custom integrators or depletion solvers can be developed by subclassing from
the following abstract base classes:

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myintegrator.rst

   abc.Integrator
   abc.SIIntegrator
   abc.DepSystemSolver

D1S Functions
-------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   d1s.prepare_tallies
   d1s.time_correction_factors
   d1s.apply_time_correction
