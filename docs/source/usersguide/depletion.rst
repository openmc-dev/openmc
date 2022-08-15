.. _usersguide_depletion:

===========================
Depletion and Transmutation
===========================

OpenMC supports transport-coupled and transport-independent depletion, or
burnup, calculations through the :mod:`openmc.deplete` Python module. OpenMC
uses transmutation reaction rates to solve a set of transmutation equations
that determine the evolution of nuclide densities within a material. The
nuclide densities predicted at some future time are then used to determine
updated reaction rates, and the process is repeated for as many timesteps as
are requested.

The depletion module is designed such that the reaction rate solution (the
transport "operator") is completely isolated from the solution of the
transmutation equations and the method used for advancing time. 

:mod:`openmc.deplete` supports multiple time-integration methods for determining
material compositions over time. Each method appears as a different class.
For example, :class:`openmc.deplete.CECMIntegrator` runs a depletion calculation
using the CE/CM algorithm (deplete over a timestep using the middle-of-step
reaction rates). An instance of :class:`~openmc.deplete.abc.TransportOperator`
is passed to one of these Integrator classes along with the timesteps and power
level::

    power = 1200.0e6  # watts
    timesteps = [10.0, 10.0, 10.0]  # days
    openmc.deplete.CECMIntegrator(op, timesteps, power, timestep_units='d').integrate()

The depletion problem is executed, and once it is done a
``depletion_results.h5`` file is written. The results can be analyzed using the
:class:`openmc.deplete.Results` class. This class has methods that allow for
easy retrieval of k-effective, nuclide concentrations, and reaction rates over
time::

    results = openmc.deplete.Results("depletion_results.h5")
    time, keff = results.get_keff()

Note that the coupling between the reaction rate solver and the transmutation
solver happens in-memory rather than by reading/writing files on disk. OpenMC
has two categories of transport operators for obtaining transmutation reaction
rates. 

.. _coupled-depletion:

Transport-coupled depletion
===========================

This category of operator solves the transport equation to obtain transmutation
reaction rates. At present, the :mod:`openmc.deplete` module offers a single
transport-coupled operator, :class:`openmc.deplete.CoupledOperator` (which uses
the OpenMC transport solver), but in principle additional transport-coupled
operator classes based on other transport codes could be implemented and no
changes to the depletion solver itself would be needed. The
:class:`openmc.deplete.CoupledOperator` class requires a :class:`~openmc.Model`
instance containing material, geometry, and settings information::

    model = openmc.Model()
    ...

    op = openmc.deplete.CoupledOperator(model)

Any material that contains a fissionable nuclide is depleted by default, but
this can behavior can be changed with the :attr:`Material.depletable` attribute.

.. important::
   
   The volume must be specified for each material that is depleted by setting
   the :attr:`Material.volume` attribute. This is necessary in order to
   calculate the proper normalization of tally results based on the source rate.

Fixed-Source Transmutation
--------------------------

When the ``power`` or ``power_density`` argument is used for one of the
Integrator classes, it is assumed that OpenMC is running in k-eigenvalue mode,
and normalization of tally results is performed based on energy deposition. It
is also possible to run a fixed-source simulation and perform normalization
based on a known source rate. First, as with all fixed-source calculations, we
need to set the run mode::

    settings.run_mode = 'fixed source'

Additionally, all materials that you wish to deplete need to be marked as such
using the :attr:`Material.depletable` attribute::

    mat = openmc.Material()
    mat.depletable = True

When constructing the :class:`~openmc.deplete.CoupledOperator`, you should
indicate that normalization of tally results will be done based on the source
rate rather than a power or power density::

    op = openmc.deplete.CoupledOperator(model, normalization_mode='source-rate')

Finally, when creating a depletion integrator, use the ``source_rates`` argument::

    integrator = openmc.deplete.PredictorIntegrator(op, timesteps, sources_rates=...)

As with the ``power`` argument, you can provide a different source rate for each
timestep in the calculation. A zero source rate for a given timestep will result
in a decay-only step, where all reaction rates are zero.

Caveats
-------

.. _energy-deposition:

Energy Deposition
~~~~~~~~~~~~~~~~~

The default energy deposition mode, ``"fission-q"``, instructs the
:class:`~openmc.deplete.CoupledOperator` to normalize reaction rates using the
product of fission reaction rates and fission Q values taken from the depletion
chain. This approach does not consider indirect contributions to energy
deposition, such as neutron heating and energy from secondary photons. In doing
this, the energy deposited during a transport calculation will be lower than
expected. This causes the reaction rates to be over-adjusted to hit the
user-specific power, or power density, leading to an over-depletion of burnable
materials.

There are some remedies. First, the fission Q values can be directly set in a
variety of ways. This requires knowing what the total fission energy release
should be, including indirect components. Some examples are provided below::

    # use a dictionary of fission_q values
    fission_q = {"U235": 202e+6}  # energy in eV

    # create a Model object
    model = openmc.Model(geometry, settings)

    # create a modified chain and write it to a new file
    chain = openmc.deplete.Chain.from_xml("chain.xml", fission_q)
    chain.export_to_xml("chain_mod_q.xml")
    op = openmc.deplete.CoupledOperator(model, "chain_mod_q.xml")

    # alternatively, pass the modified fission Q directly to the operator
    op = openmc.deplete.CoupledOperator(model, "chain.xml",
        fission_q=fission_q)


A more complete way to model the energy deposition is to use the modified
heating reactions described in :ref:`methods_heating`. These values can be used
to normalize reaction rates instead of using the fission reaction rates with::

    op = openmc.deplete.CoupledOperator(model, "chain.xml",
        normalization_mode="energy-deposition")

These modified heating libraries can be generated by running the latest version
of :meth:`openmc.data.IncidentNeutron.from_njoy()`, and will eventually be bundled
into the distributed libraries.

Local Spectra and Repeated Materials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is not uncommon to explicitly create a single burnable material across many
locations. From a pure transport perspective, there is nothing wrong with
creating a single 3.5 wt.% enriched fuel ``fuel_3``, and placing that fuel in
every fuel pin in an assembly or even full core problem. This certainly
expedites the model making process, but can pose issues with depletion. Under
this setup, :mod:`openmc.deplete` will deplete a single ``fuel_3`` material
using a single set of reaction rates, and produce a single new composition for
the next time step. This can be problematic if the same ``fuel_3`` is used in
very different regions of the problem.

As an example, consider a full-scale power reactor core with vacuum boundary
conditions, and with fuel pins solely composed of the same ``fuel_3`` material.
The fuel pins towards the center of the problem will surely experience a more
intense neutron flux and greater reaction rates than those towards the edge of
the domain. This indicates that the fuel in the center should be at a more
depleted state than periphery pins, at least for the fist depletion step.
However, without any other instructions, OpenMC will deplete ``fuel_3`` as a
single material, and all of the fuel pins will have an identical composition at
the next transport step.

This can be countered by instructing the operator to treat repeated instances
of the same material as a unique material definition with::

    op = openmc.deplete.CoupledOperator(model, chain_file,
        diff_burnable_mats=True)

For our example problem, this would deplete fuel on the outer region of the
problem with different reaction rates than those in the center. Materials will
be depleted corresponding to their local neutron spectra, and have unique
compositions at each transport step.  The volume of the original ``fuel_3``
material must represent the volume of **all** the ``fuel_3`` in the problem.
When creating the unique materials, this volume will be equally distributed
across all material instances.


.. note::

    This will increase the total memory usage and run time due to an increased
    number of tallies and material definitions.

Transport-independent depletion
===============================

.. warning::
   
   This feature is still under heavy development and has yet to be rigorously 
   verified. API changes and feature additions are possible and likely in
   the near future.

This category of operator uses pre-calculated one-group microscopic cross
sections to obtain transmutation reaction rates. OpenMC provides the
:class:`~openmc.deplete.IndependentOperator` for this method of calculation.
While the one-group microscopic cross sections can be calculated using a
transport solver, :class:`~openmc.deplete.IndependentOperator` is not directly
coupled to any transport solver. The
:class:`~openmc.deplete.IndependentOperator` class requires a
:class:`openmc.Materials` object, a :class:`~openmc.deplete.MicroXS` object,
and a path to a depletion chain file::

    # load in the microscopic cross sections
    materials = openmc.Materials()
    ...

    micro_xs = openmc.deplete.MicroXS.from_csv(micro_xs_path)
    op = openmc.deplete.IndependentOperator(materials, micro_xs, chain_file)

.. note::

   The same statements from :ref:`coupled-depletion` about which
   materials are depleted and the requirement for depletable materials to have
   a specified volume also apply here.

An alternate constructor,
:meth:`~openmc.deplete.IndependentOperator.from_nuclides`, accepts a volume and
dictionary of nuclide concentrations in place of the :class:`openmc.Materials`
object::

    nuclides = {'U234': 8.92e18,
                'U235': 9.98e20,
                'U238': 2.22e22,
                'U236': 4.57e18,
                'O16': 4.64e22,
                'O17': 1.76e19}
    volume = 0.5 
    op = openmc.deplete.IndependentOperator.from_nuclides(volume,
                                                          nuclides,
                                                          micro_xs,
                                                          chain_file,
                                                          nuc_units='atom/cm3')

A user can then define an integrator class as they would for a coupled
transport-depletion calculation and follow the same steps from there.

.. note::

   Ideally, one-group cross section data should be available for every
   reaction in the depletion chain. If a nuclide that has a reaction 
   associated with it in the depletion chain is present in the `nuclides` 
   parameter but not the cross section data, that reaction will not be
   simulated.

Generating Microscopic Cross Sections
-------------------------------------

Users can generate the one-group microscopic cross sections needed by
:class:`~openmc.deplete.IndependentOperator` using the
:class:`~openmc.deplete.MicroXS` class::

    import openmc

    model = openmc.Model.from_xml()

    micro_xs = openmc.deplete.MicroXS.from_model(model,
                                                 model.materials[0],
                                                 chain_file)

The :meth:`~openmc.deplete.MicroXS.from_model()` method will produce a
:class:`~openmc.deplete.MicroXS` object with microscopic cross section data in
units of barns, which is what :class:`~openmc.deplete.IndependentOperator`
expects the units to be. The :class:`~openmc.deplete.MicroXS` class also
includes functions to read in cross section data directly from a ``.csv`` file
or from data arrays::

    micro_xs = MicroXS.from_csv(micro_xs_path)

    nuclides = ['U234', 'U235', 'U238']
    reactions = ['fission', '(n,gamma)']
    data = np.array([[0.1, 0.2],
                     [0.3, 0.4],
                     [0.01, 0.5]])
    micro_xs = MicroXS.from_array(nuclides, reactions, data)

.. important::

   Both :meth:`~openmc.deplete.MicroXS.from_csv()` and
   :meth:`~openmc.deplete.MicroXS.from_array()` assume the cross section values 
   provided are in barns by defualt, but have no way of verifying this. Make
   sure your cross sections are in the correct units before passing to a
   :class:`~openmc.deplete.IndependentOperator` object.

Caveats
-------

Reaction Rate Normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`~openmc.deplete.IndependentOperator` class supports two methods for
normalizing reaction rates:

.. important::

   Make sure you set the correct parameter in the :class:`openmc.abc.Integrator`
   class. Use the ``source_rates`` parameter when
   ``normalization_mode == source-rate``, and use ``power`` or ``power_density``
   when ``normalization_mode == fission-q``.

1. ``source-rate`` normalization, which assumes the ``source_rate`` provided by
   the time integrator is a flux, and obtains the reaction rates by multiplying
   the cross-sections by the ``source-rate``.
2. ``fission-q`` normalization, which uses the ``power`` or ``power_density``
   provided by the time integrator to obtain reaction rates by computing a value
   for the flux based on this power. The general equation for the flux is 

   .. math::

      \phi = \frac{P}{V \cdot \sum_i (Q_i \cdot \sigma^f_i \cdot \n_i)}

   where :math:`\sum_i` is the sum over all nuclides :math:`i`. This equation
   makes the same assumptions and issues as discussed in
   :ref:`energy-deposition`. Unfortunately, the proposed solution in that
   section does not apply here since we are decoupled from transport code.
   However, there is a method to converge to a more accurate value for flux by
   using substeps during time integration.
   `This paper <https://doi.org/10.1016/j.anucene.2016.05.031>`_ provides a
   good discussion of this method. 

.. warning::

   The accuracy of results when using ``fission-q`` is entirely dependent on
   your depletion chain. Make sure it has sufficient data to resolve the
   dynamics of your particular scenario. 

Multiple Materials
~~~~~~~~~~~~~~~~~~

Running a depletion simulation with multiple materials using the
``source-rate`` normalization method treats each material as completely
separate with respect to reaction rates. This can be useful for running many
different cases of a particular scenario. However, running a depletion
simulation with multiple materials using the ``fission-q`` normalization method
treats each material as part of the same "reactor" due to how ``fission-q``
normalization accumulates energy values from each material to a single value.
This behavior may change in the future.

Time integration
~~~~~~~~~~~~~~~~

The one-group microscopic cross sections passed to
:class:`openmc.deplete.IndependentOperator` are fixed values for the entire
depletion simulation. This implicit assumption may produce inaccurate results
for certain scenarios.
