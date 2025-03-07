import matplotlib.pyplot as plt
import numpy as np
import openmc


# Get results from statepoint
with openmc.StatePoint("statepoint.10.h5") as sp:
    tally = sp.get_tally(name="Density")

# Get the time grid
t = tally.filters[0].values
dt = np.diff(t)
t_mid = 0.5 * (t[:-1] + t[1:])

# Bin-averaged result
density_mean = tally.mean.ravel() / dt

# Plot flux spectrum
fig, ax = plt.subplots()
ax.loglog(t_mid, density_mean, "ok", fillstyle="none")
ax.set_xlabel("Time [s]")
ax.set_ylabel("Total density")
ax.grid()
plt.show()
