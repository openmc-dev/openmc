import openmc.examples
import numpy as np
import matplotlib.pyplot as plt

model = openmc.examples.pwr_pin_cell()

model.tallies

# Create equal-lethargy energies to put in filter
energies = np.logspace(np.log10(1e-5), np.log10(20.0e6), 501)
e_filter = openmc.EnergyFilter(energies)

# Create tally with energy filter
tally = openmc.Tally()
tally.filters = [e_filter]
tally.scores = ['flux']

# Set model tallies
model.tallies = [tally]

openmc.mgxs.GROUP_STRUCTURES.keys()

# Create energy filter using SHEM-361 group structure
energies_shem = openmc.mgxs.GROUP_STRUCTURES['SHEM-361']
shem_filter = openmc.EnergyFilter(openmc.mgxs.GROUP_STRUCTURES['SHEM-361'])

tally_shem = openmc.Tally()
tally_shem.filters = [shem_filter]
tally_shem.scores = ['flux']

model.tallies.append(tally_shem)

model.settings.particles = 10000
model.settings.batches = 50

model.settings.statepoint = {"overwrite_latest" : "2"}

print(model.settings.statepoint)
sp_path = model.run(output=True)

with openmc.StatePoint(sp_path) as sp:
    t = sp.tallies[tally.id]
    flux500_mean = t.mean.ravel()
    flux500_unc = t.std_dev.ravel()
    
    t_shem = sp.tallies[tally_shem.id]
    flux_shem_mean = t_shem.mean.ravel()
    flux_shem_unc = t_shem.std_dev.ravel()

    fig, ax = plt.subplots()

ax.step(energies[:-1], flux500_mean/np.diff(energies), where='post', label='500 group')
ax.step(energies_shem[:-1], flux_shem_mean/np.diff(energies_shem), where='post', label='SHEM-361')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('Flux [n-cm/eV-src]')
ax.grid()
ax.legend()


fig, ax = plt.subplots()
ax.loglog(energies[:-1], flux500_mean, '.', color='C0', label='500 group')
ax.loglog(energies_shem[:-1], flux_shem_mean, '.', color='C1', label='SHEM-361')
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('Flux [n-cm/src]')
ax.grid()
ax.legend()