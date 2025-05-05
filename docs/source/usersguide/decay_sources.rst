.. usersguide_decay_sources:

=============
Decay Sources
=============

Through the :ref:`depletion <usersguide_depletion>` capabilities in OpenMC, it
is possible to simulate radiation emitted from the decay of activated materials.
For fusion energy systems, this is commonly done using what is known as the
`rigorous 2-step <https://doi.org/10.1016/S0920-3796(02)00144-8>`_ (R2S) method.
In this method, a neutron transport calculation is used to determine the neutron
flux and reaction rates over a cell- or mesh-based spatial discretization of the
model. Then, the neutron flux in each discrete region is used to predict the
activated material composition using a depletion solver. Finally, a photon
transport calculation with a source based on the activity and energy spectrum of
the activated materials is used to determine a desired physical response (e.g.,
a dose rate) at one or more locations of interest.

Once a depletion simulation has been completed in OpenMC, the intrinsic decay
source can be determined as follows. First the activated material composition
can be determined using the :class:`openmc.deplete.Results` object. Indexing an
instance of this class with the timestep index returns a
:class:`~openmc.deplete.StepResult` object, which itself has a
:meth:`~openmc.deplete.StepResult.get_material` method. Once the activated
:class:`~openmc.Material` has been obtained, the
:meth:`~openmc.Material.get_decay_photon_energy` method will give the energy
spectrum of the decay photon source. The integral of the spectrum also indicates
the intensity of the source in units of [Bq]. Altogether, the workflow looks as
follows::

    results = openmc.deplete.Results("depletion_results.h5")

    # Get results at last timestep
    step = results[-1]

    # Get activated material composition for ID=1
    activated_mat = step.get_material('1')

    # Determine photon source
    photon_energy = activated_mat.get_decay_photon_energy()

By default, the :meth:`~openmc.Material.get_decay_photon_energy` method will
eliminate spectral lines with very low intensity, but this behavior can be
configured with the ``clip_tolerance`` argument.

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

