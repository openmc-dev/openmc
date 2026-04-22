.. _random_ray:

===================
Kinetic Simulations
===================

Much of OpenMC's existing infrastructure assumes a system is at a steady state.
Users interested in reactor transients should use the kinetic simulation
capability. Kinetic simulations can be enabled as::

    settings.kinetic_simulation = True

Kinetic simulations require the user to specify time step settings using the
``settings.timestep_parameters`` object. As an example, a 10 second long
simulation using 1 millisecond time steps can be set as::

    settings.timestep_parameters = {'n_timesteps': 1000,
                                    'dt': 1,
                                    'timestep_units': 'ms'}

.. note::
    When running kineic simulations, OpenMC will generate a statepoint file for
    each time step named ``openmc_td_simulation_{i}.h5``, where ``i`` is the
    index of the time step.

Currently, only material density transients can be simulated. A density
transient can be specified using the Python API::

   water = openmc.Material(name='water')
    
   # Linear ramp transient
   densities = np.linspace(1, 0.95, 1000)

   water.set_density('macro', 1.0, density_timeseries=densities)

.. note::
    Kinetic simulations are currently only supported in eigenvalue mode using the Random Ray solver.

----------------------
Random Ray Quick Start
----------------------

Random Ray's :ref:`automatic setup workflow <quick_start>` can
be utilized to quickly generate models for kinetic simulations, however
slight modifications are needed:

1. The `kinetic` flag should be set to `True`, and the number of delay groups must be
   specified via `num_delay_groups` in the call to :meth:`openmc.Model.convert_to_multigroup`.
   Take care to not specify more delay groups than the cross section library you
   are using supports, as this will cause all sorts of problems when you try to run
   the simulation.
2. By default, a time step size of 1 ms is chosen to accurately resolve delayed
   neutron precursor dynamics, but the number of time steps must be manually
   set in `openmc.Settings.timestep_parameters['n_timesteps']`.
3. Material density timeseries' must be manually added.

An example process of converting an existing continuous energy
Monte Carlo model to a random ray model for kinetic simulations
is show below::

  # Define continuous energy model as normal
  model = openmc.Model()
  ...

  # Convert model to kinetic multigroup (will auto-generate MGXS library if needed)
  # Most cross-section libraries support 6 delay groups
  model.convert_to_multigroup(kinetic=True, num_delayed_groups=6)

  # Add required manual kinetic simulation settigns
  model.settings.timestep_parameters['n_timesteps'] = 5

  # Add a material density timeseries to material 0 in the model
  density_timeseries = np.linspace(1, 0.95, 100)
  model.materials[0].set_density('macro', density=1.0, density_timeseries=density_timeseries)

  # Convert model to random ray and initialize random ray parameters
  # to reasonable defaults based on the specifics of the geometry
  model.convert_to_random_ray()

  # (Optional) Overlay source region decomposition mesh to improve fidelity of the
  # random ray solver. Adjust 'n' for fidelity vs runtime.
  n = 100
  mesh = openmc.RegularMesh()
  mesh.dimension = (n, n, n)
  mesh.lower_left = model.geometry.bounding_box.lower_left
  mesh.upper_right = model.geometry.bounding_box.upper_right
  model.settings.random_ray['source_region_meshes'] = [(mesh, [model.geometry.root_universe])]

  # (Optional) Increase the number of rays/batch, to reduce uncertainty
  model.settings.particles = 500
