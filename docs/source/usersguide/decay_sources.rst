.. usersguide_decay_sources:

=============
Decay Sources
=============

Through the :ref:`depletion <usersguide_depletion>` capabilities in OpenMC, it
is possible to simulate radiation emitted from the decay of activated materials.
For fusion energy systems, this is commonly done using either the `rigorous
2-step <https://doi.org/10.1016/S0920-3796(02)00144-8>`_ (R2S) method or the
`direct 1-step <https://doi.org/10.1016/S0920-3796(01)00188-0>`_ (D1S) method.
In the R2S method, a neutron transport calculation is used to determine the
neutron flux and reaction rates over a cell- or mesh-based spatial
discretization of the model. Then, the neutron flux in each discrete region is
used to predict the activated material composition using a depletion solver.
Finally, a photon transport calculation with a source based on the activity and
energy spectrum of the activated materials is used to determine a desired
physical response (e.g., a dose rate) at one or more locations of interest.
OpenMC includes automation for both the R2S and D1S methods as described in the
following sections.

Rigorous 2-Step (R2S) Calculations
==================================

OpenMC includes an :class:`openmc.deplete.R2SManager` class that fully automates
cell- and mesh-based R2S calculations. Before we describe this class, it is
useful to understand the basic mechanics of how an R2S calculation works.
Generally, it involves the following steps:

1. The :meth:`openmc.deplete.get_microxs_and_flux` function is called to run a
   neutron transport calculation that determines fluxes and microscopic cross
   sections in each activation region.
2. The :class:`openmc.deplete.IndependentOperator` and
   :class:`openmc.deplete.PredictorIntegrator` classes are used to carry out a
   depletion (activation) calculation in order to determine predicted material
   compositions based on a set of timesteps and source rates.
3. The activated material composition is determined using the
   :class:`openmc.deplete.Results` class. Indexing an instance of this class
   with the timestep index returns a :class:`~openmc.deplete.StepResult` object,
   which itself has a :meth:`~openmc.deplete.StepResult.get_material` method
   returning an activated material.
4. The :meth:`openmc.Material.get_decay_photon_energy` method is used to obtain
   the energy spectrum of the decay photon source. The integral of the spectrum
   also indicates the intensity of the source in units of [Bq].
5. A new photon source is defined using one of OpenMC's source classes with the
   energy distribution set equal to the object returned by the
   :meth:`openmc.Material.get_decay_photon_energy` method. The source is then
   assigned to a photon :class:`~openmc.Model`.
6. A photon transport calculation is run with ``model.run()``.

Altogether, the workflow looks as follows::

    # Run neutron transport calculation
    fluxes, micros = openmc.deplete.get_microxs_and_flux(model, domains)

    # Run activation calculation
    op = openmc.deplete.IndependentOperator(mats, fluxes, micros)
    timesteps = ...
    source_rates = ...
    integrator = openmc.deplete.Integrator(op, timesteps, source_rates)
    integrator.integrate()

    # Get decay photon source at last timestep
    results = openmc.deplete.Results("depletion_results.h5")
    step = results[-1]
    activated_mat = step.get_material('1')
    photon_energy = activated_mat.get_decay_photon_energy()
    photon_source = openmc.IndependentSource(
        space=...,
        energy=photon_energy,
        particle='photon',
        strength=photon_energy.integral()
    )

    # Run photon transport calculation
    model.settings.source = photon_source
    model.run()

Note that by default, the :meth:`~openmc.Material.get_decay_photon_energy`
method will eliminate spectral lines with very low intensity, but this behavior
can be configured with the ``clip_tolerance`` argument.

Cell-based R2S
--------------

In practice, users do not need to manually go through each of the steps in an R2S
calculation described above. The :class:`~openmc.deplete.R2SManager` fully
automates the execution of neutron transport, depletion, decay source
generation, and photon transport. For a cell-based R2S calculation, once you
have a :class:`~openmc.Model` that has been defined, simply create an instance
of :class:`~openmc.deplete.R2SManager` by passing the model and a list of cells
to activate::

    r2s = openmc.deplete.R2SManager(model, [cell1, cell2, cell3])

Note that the ``volume`` attribute must be set for any cell that is to be
activated. The :class:`~openmc.deplete.R2SManager` class allows you to
optionally specify a separate photon model; if not given as an argument, it will
create a shallow copy of the original neutron model (available as the
``neutron_model`` attribute) and store it in the ``photon_model`` attribute. We
can use this to define tallies specific to the photon model::

    dose_tally = openmc.Tally()
    ...
    r2s.photon_model.tallies = [dose_tally]

Next, define the timesteps and source rates for the activation calculation::

    timesteps = [(3.0, 'd'), (5.0, 'h')]
    source_rates = [1e12, 0.0]

In this case, the model is irradiated for 3 days with a source rate of
:math:`10^{12}` neutron/sec and then the source is turned off and the activated
materials are allowed to decay for 5 hours. These parameters should be passed to
the :meth:`~openmc.deplete.R2SManager.run` method to execute the full R2S
calculation. Before we can do that though, for a cell-based calculation, the one
other piece of information that is needed is bounding boxes of the activated
cells::

    bounding_boxes = {
        cell1.id: cell1.bounding_box,
        cell2.id: cell2.bounding_box,
        cell3.id: cell3.bounding_box
    }

Note that calling the ``bounding_box`` attribute may not work for all
constructive solid geometry regions (for example, a cell that uses a
non-axis-aligned plane). In these cases, the bounding box will need to be
specified manually. Once you have a set of bounding boxes, the R2S calculation
can be run::

    r2s.run(timesteps, source_rates, bounding_boxes=bounding_boxes)

If not specified otherwise, a photon transport calculation is run at each time
in the depletion schedule. That means in the case above, we would see three
photon transport calculations. To specify specific times at which photon
transport calculations should be run, pass the ``photon_time_indices`` argument.
For example, if we wanted to run a photon transport calculation only on the last
time (after the 5 hour decay), we would run::

    r2s.run(timesteps, source_rates, bounding_boxes=bounding_boxes,
            photon_time_indices=[2])

After an R2S calculation has been run, the :class:`~openmc.deplete.R2SManager`
instance will have a ``results`` dictionary that allows you to directly access
results from each of the steps. It will also write out all the output files into
a directory that is named "r2s_<timestamp>/". The ``output_dir`` argument to the
:meth:`~openmc.deplete.R2SManager.run` method enables you to override the
default output directory name if desired.

The :meth:`~openmc.deplete.R2SManager.run` method actually runs three
lower-level methods under the hood::

    r2s.step1_neutron_transport(...)
    r2s.step2_activation(...)
    r2s.step3_photon_transport(...)

For users looking for more control over the calculation, these lower-level
methods can be used in lieu of the :meth:`openmc.deplete.R2SManager.run` method.

Mesh-based R2S
--------------

Executing a mesh-based R2S calculation looks nearly identical to the cell-based
R2S workflow described above. The only difference is that instead of passing a
list of cells to the ``domains`` argument of
:class:`~openmc.deplete.R2SManager`, you need to define a mesh object and pass
that instead. This might look like the following::

    # Define a regular Cartesian mesh
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-50., -50., 0.)
    mesh.upper_right = (50., 50., 75.)
    mesh.dimension = (10, 10, 5)

    r2s = openmc.deplete.R2SManager(model, mesh)

Executing the R2S calculation is then performed by adding photon tallies and
calling the :meth:`~openmc.deplete.R2SManager.run` method with the appropriate
timesteps and source rates. Note that in this case we do not need to define cell
volumes or bounding boxes as is required for a cell-based R2S calculation.
Instead, during the neutron transport step, OpenMC will run a raytracing
calculation to determine material volume fractions within each mesh element
using the :meth:`openmc.MeshBase.material_volumes` method. Arguments to this
method can be customized via the ``mat_vol_kwargs`` argument to the
:meth:`~openmc.deplete.R2SManager.run` method. Most often, this would involve
customizing the number of rays traced to obtain better estimates of volumes. As
an example, if we wanted to run the raytracing calculation with 10 million rays,
we would run::

    r2s.run(timesteps, source_rates, mat_vol_kwargs={'n_samples': 10_000_000})

Direct 1-Step (D1S) Calculations
================================

OpenMC also includes built-in capability for performing shutdown dose rate
calculations using the `direct 1-step
<https://doi.org/10.1016/S0920-3796(01)00188-0>`_ (D1S) method. In this method,
a single coupled neutron--photon transport calculation is used where the prompt
photon production is replaced with photons produced from the decay of
radionuclides in an activated material. To obtain properly scaled results, it is
also necessary to apply time correction factors. A normal neutron transport
calculation can be extended to a D1S calculation with a few helper functions.
First, import the ``d1s`` submodule, which is part of :mod:`openmc.deplete`::

    from openmc.deplete import d1s

First, you need to instruct OpenMC to use decay photon data instead of prompt
photon data. This is done with an attribute on the :class:`~openmc.Settings`
class::

    model = openmc.Model()
    ...
    model.settings.use_decay_photons = True

To prepare any tallies for use of the D1S method, you should call the
:func:`~openmc.deplete.d1s.prepare_tallies` function, which adds a
:class:`openmc.ParentNuclideFilter` (used later for assigning time correction
factors) to any applicable tally and returns a list of possible radionuclides
based on the :ref:`chain file <usersguide_data>`. Once the tallies are prepared,
the model can be simulated::

    output_path = model.run()

Finally, the time correction factors need to be computed and applied to the
relevant tallies. This can be done with the aid of the
:func:`~openmc.deplete.d1s.time_correction_factors` and
:func:`~openmc.deplete.d1s.apply_time_correction` functions::

    # Compute time correction factors based on irradiation schedule
    factors = d1s.time_correction_factors(nuclides, timesteps, source_rates)

    # Get tally from statepoint
    with openmc.StatePoint(output_path) as sp:
        dose_tally = sp.get_tally(name='dose tally')

    # Apply time correction factors
    tally = d1s.apply_time_correction(dose_tally, factors, time_index)

