.. _usersguide_depletion:

=========
Depletion
=========

OpenMC supports coupled depletion, or burnup, calculations through the
:mod:`openmc.deplete` Python module. OpenMC solves the transport equation to
obtain transmutation reaction rates, and then the reaction rates are used to
solve a set of transmutation equations that determine the evolution of nuclide
densities within a material. The nuclide densities predicted as some future time
are then used to determine updated reaction rates, and the process is repeated
for as many timesteps as are requested.

The depletion module is designed such that the flux/reaction rate solution (the
transport "operator") is completely isolated from the solution of the
transmutation equations and the method used for advancing time. At present, the
:mod:`openmc.deplete` module offers a single transport operator,
:class:`openmc.deplete.Operator` (which uses the OpenMC transport solver), but
in principle additional operator classes based on other transport codes could be
implemented and no changes to the depletion solver itself would be needed. The
operator class requires a :class:`openmc.Geometry` instance and a
:class:`openmc.Settings` instance::

    geom = openmc.Geometry()
    settings = openmc.Settings()
    ...

    op = openmc.deplete.Operator(geom, settings)

:mod:`openmc.deplete` supports multiple time-integration methods for determining
material compositions over time. Each method appears as a different function.
For example, :func:`openmc.deplete.integrator.cecm` runs a depletion calculation
using the CE/CM algorithm (deplete over a timestep using the middle-of-step
reaction rates). An instance of :class:`openmc.deplete.Operator` is passed to
one of these functions along with the power level and timesteps::

    power = 1200.0e6
    days = 24*60*60
    timesteps = [10.0*days, 10.0*days, 10.0*days]
    openmc.deplete.cecm(op, power, timesteps)

The coupled transport-depletion problem is executed, and once it is done a
``depletion_results.h5`` file is written. The results can be analyzed using the
:class:`openmc.deplete.ResultsList` class. This class has methods that allow for
easy retrieval of k-effective, nuclide concentrations, and reaction rates over
time::

    results = openmc.deplete.ResultsList("depletion_results.h5")
    time, keff = results.get_eigenvalue()

Note that the coupling between the transport solver and the transmutation solver
happens in-memory rather than by reading/writing files on disk.
