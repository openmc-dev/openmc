import matplotlib.pyplot as plt
import openmc

# Get the flux from the statepoint
with openmc.StatePoint('statepoint.10.h5') as sp:
    flux = sp.tallies[1].mean
    flux.shape = (50, 50)

# Plot the flux
fig, ax = plt.subplots()
ax.imshow(flux, origin='lower', extent=(-5.0, 5.0, -5.0, 5.0))
ax.set_xlabel('x [cm]')
ax.set_ylabel('y [cm]')
plt.show()

# If all worked well, you should see a ring "imprint" as well as a higher flux
# to the right side (since the custom source has all particles moving in the
# positive x direction)
