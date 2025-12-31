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
