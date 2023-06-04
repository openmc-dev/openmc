import matplotlib.pyplot as plt
import openmc
import numpy
import os

averages = numpy.ndarray(shape=(0))
particles = numpy.ndarray(shape=(0))

for nparticles in range(100, 1000, 100):
    os.system("python sim.py " + str(nparticles))  

    # Get results from statepoint
    with openmc.StatePoint('statepoint.100.h5') as sp:
        t = sp.get_tally(name="Flux spectrum")

        # Get the energies from the energy filter
        energy_filter = t.filters[0]
        energies = energy_filter.bins[:, 0]

        # Get the flux values
        mean = t.get_values(value='mean').ravel()
        uncertainty = t.get_values(value='std_dev').ravel()

    averages = numpy.append(averages, numpy.average(mean))
    particles = numpy.append(particles, nparticles)

# Plot flux spectrum
fix, ax = plt.subplots()
ax.loglog(particles, averages, drawstyle='steps-post')
ax.set_xlabel('Number of Particles Sent Into settings.particles')
ax.set_ylabel('Average Flux for the Whole Run')
ax.grid(True, which='both')
plt.show()
