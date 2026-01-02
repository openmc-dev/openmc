import openmc
import openmc.deplete
import matplotlib.pyplot as plt

###############################################################################
#                      Load previous simulation results
###############################################################################

# Load geometry from statepoint
statepoint = 'statepoint.100.h5'
with openmc.StatePoint(statepoint) as sp:
    geometry = sp.summary.geometry

# Load previous depletion results
previous_results = openmc.deplete.Results("depletion_results.h5")

###############################################################################
#                      Transport calculation settings
###############################################################################

# Instantiate a Settings object, set all runtime parameters
settings = openmc.Settings()
settings.batches = 100
settings.inactive = 10
settings.particles = 10000

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-0.62992, -0.62992, -1, 0.62992, 0.62992, 1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
settings.source = openmc.IndependentSource(
    space=uniform_dist, constraints={'fissionable': True})

entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left = [-0.39218, -0.39218, -1.e50]
entropy_mesh.upper_right = [0.39218, 0.39218, 1.e50]
entropy_mesh.dimension = [10, 10, 1]
settings.entropy_mesh = entropy_mesh

###############################################################################
#                   Initialize and run depletion calculation
###############################################################################

model = openmc.Model(geometry=geometry, settings=settings)

# Create depletion "operator"
chain_file = 'chain_simple.xml'
op = openmc.deplete.CoupledOperator(model, chain_file, previous_results)

# Perform simulation using the predictor algorithm
time_steps = [1.0, 1.0, 1.0, 1.0, 1.0]  # days
power = 174  # W/cm, for 2D simulations only (use W for 3D)
integrator = openmc.deplete.PredictorIntegrator(op, time_steps, power, timestep_units='d')
integrator.integrate()

###############################################################################
#                    Read depletion calculation results
###############################################################################

# Open results file
results = openmc.deplete.Results("depletion_results.h5")

# Obtain K_eff as a function of time
time, keff = results.get_keff(time_units='d')

# Obtain U235 concentration as a function of time
uo2 = geometry.get_all_material_cells()[1]
_, n_U235 = results.get_atoms(uo2, 'U235')

# Obtain Xe135 capture reaction rate as a function of time
_, Xe_capture = results.get_reaction_rate(uo2, 'Xe135', '(n,gamma)')

###############################################################################
#                            Generate plots
###############################################################################

fig, ax = plt.subplots()
ax.errorbar(time, keff[:, 0], keff[:, 1], label="K-effective")
ax.set_xlabel("Time [d]")
ax.set_ylabel("Keff")
plt.show()

fig, ax = plt.subplots()
ax.plot(time, n_U235, label="U235")
ax.set_xlabel("Time [d]")
ax.set_ylabel("U235 atoms")
plt.show()

fig, ax = plt.subplots()
ax.plot(time, Xe_capture, label="Xe135 capture")
ax.set_xlabel("Time [d]")
ax.set_ylabel("Xe135 capture rate")
plt.show()
