#!/usr/bin/env python3
"""Plot FNS flux spectrum: CCFE-709 vs converted UKAEA-1102."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
from pathlib import Path
import sys

import openmc.mgxs

# Load FNS 709-group flux
flux_file = 'fns_flux_709.npy'
flux_709  = np.load(flux_file)

# Get energy structures
E_709  = openmc.mgxs.GROUP_STRUCTURES['CCFE-709']
E_1102 = openmc.mgxs.GROUP_STRUCTURES['UKAEA-1102']

# Convert flux
flux_1102 = openmc.mgxs.convert_flux_groups(flux_709, 'CCFE-709', 'UKAEA-1102')

# Calculate flux per unit lethargy: phi/u = phi / ln(E_hi/E_lo)
def flux_per_lethargy(flux, edges):
    lethargy = np.log(edges[1:] / edges[:-1])
    return flux / lethargy

flux_709_per_u = flux_per_lethargy(flux_709, E_709)
flux_1102_per_u = flux_per_lethargy(flux_1102, E_1102)

# Energy range
E_min, E_max = 1e-3, 20e6

# Create figure with 2x2 layout (left 67%, right 33%)
fig, axes = plt.subplots(2, 2, figsize=(16, 10), gridspec_kw={'width_ratios': [2, 1]})
ax1, ax_thermal = axes[0]
ax2, ax_fast = axes[1]

# Top-left: Flux vs Energy (full range)
ax1.stairs(flux_709, E_709, label=f'CCFE-709 ({len(flux_709)} groups)', color='blue', linewidth=1.2)
ax1.stairs(flux_1102, E_1102, label=f'UKAEA-1102 ({len(flux_1102)} groups)', color='red', linewidth=0.8, alpha=0.8)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(E_min, E_max)
ax1.set_ylabel('Flux [n/cm²/s]')
ax1.set_title('FNS Neutron Flux: CCFE-709 vs UKAEA-1102 Conversion')
ax1.legend(loc='upper left')
ax1.grid(True, alpha=0.3, which='both')
ax1.xaxis.set_major_locator(LogLocator(base=10, numticks=15))
ax1.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))

# Bottom-left: Flux per unit lethargy vs Energy (full range)
ax2.stairs(flux_709_per_u, E_709, label=f'CCFE-709 ({len(flux_709)} groups)', color='blue', linewidth=1.2)
ax2.stairs(flux_1102_per_u, E_1102, label=f'UKAEA-1102 ({len(flux_1102)} groups)', color='red', linewidth=0.8, alpha=0.8)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(E_min, E_max)
ax2.set_xlabel('Energy [eV]')
ax2.set_ylabel('Flux per unit lethargy [n/cm²/s]')
ax2.legend(loc='upper left')
ax2.grid(True, alpha=0.3, which='both')
ax2.xaxis.set_major_locator(LogLocator(base=10, numticks=15))
ax2.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))

# Helper to get y-limits for a given energy range
def get_ylim_for_range(flux_per_u, edges, e_min, e_max):
    mask = (edges[:-1] >= e_min) & (edges[1:] <= e_max)
    vals = flux_per_u[mask]
    vals = vals[vals > 0]
    if len(vals) == 0:
        return None, None
    return vals.min() * 0.8, vals.max() * 1.2

# Top-right: Thermal zoom (1-10 eV) - flux per lethargy
thermal_emin, thermal_emax = 1, 10
ax_thermal.stairs(flux_709_per_u, E_709, label='CCFE-709', color='blue', linewidth=1.5)
ax_thermal.stairs(flux_1102_per_u, E_1102, label='UKAEA-1102', color='red', linewidth=1.0, alpha=0.8)
ax_thermal.set_xscale('log')
ax_thermal.set_yscale('log')
ax_thermal.set_xlim(thermal_emin, thermal_emax)
ymin1, ymax1 = get_ylim_for_range(flux_709_per_u, E_709, thermal_emin, thermal_emax)
ymin2, ymax2 = get_ylim_for_range(flux_1102_per_u, E_1102, thermal_emin, thermal_emax)
if ymin1 and ymin2:
    ax_thermal.set_ylim(min(ymin1, ymin2), max(ymax1, ymax2))
ax_thermal.set_ylabel('Flux/lethargy [n/cm²/s]')
ax_thermal.set_title('Thermal Zoom (1-10 eV)')
ax_thermal.legend(loc='best', fontsize=8)
ax_thermal.grid(True, alpha=0.3, which='both')

# Bottom-right: Fast zoom (6-16 MeV) - flux per lethargy
fast_emin, fast_emax = 6e6, 16e6
ax_fast.stairs(flux_709_per_u, E_709, label='CCFE-709', color='blue', linewidth=1.5)
ax_fast.stairs(flux_1102_per_u, E_1102, label='UKAEA-1102', color='red', linewidth=1.0, alpha=0.8)
ax_fast.set_xscale('log')
ax_fast.set_yscale('log')
ax_fast.set_xlim(fast_emin, fast_emax)
ymin1, ymax1 = get_ylim_for_range(flux_709_per_u, E_709, fast_emin, fast_emax)
ymin2, ymax2 = get_ylim_for_range(flux_1102_per_u, E_1102, fast_emin, fast_emax)
if ymin1 and ymin2:
    ax_fast.set_ylim(min(ymin1, ymin2), max(ymax1, ymax2))
ax_fast.set_xlabel('Energy [eV]')
ax_fast.set_ylabel('Flux/lethargy [n/cm²/s]')
ax_fast.set_title('Fast Zoom (6-16 MeV)')
ax_fast.legend(loc='best', fontsize=8)
ax_fast.grid(True, alpha=0.3, which='both')

# Add conservation check
total_709 = np.sum(flux_709)
total_1102 = np.sum(flux_1102)
ax1.text(0.98, 0.02, f'Total flux: 709={total_709:.4e}, 1102={total_1102:.4e}\nRelative diff: {abs(total_709-total_1102)/total_709:.2e}',
         transform=ax1.transAxes, fontsize=9, verticalalignment='bottom', horizontalalignment='right',
         family='monospace', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig(Path(__file__).parent / 'fns_flux_conversion_comparison.png', dpi=500)
print(f"Saved: {Path(__file__).parent / 'fns_flux_conversion_comparison.png'}")
plt.show()
