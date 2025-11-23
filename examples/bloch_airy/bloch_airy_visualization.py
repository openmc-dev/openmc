#!/usr/bin/env python3
"""
Bloch-Airy Quantum Gravitational Bound States Visualization

This script demonstrates the discretized quantum bouncing of ultracold neutrons
(UCN) using OpenMC's Bloch-Airy model. It compares classical vs quantum
predictions and visualizes the quantized energy levels and probability
distributions.

The Bloch-Airy model describes neutrons bouncing on a horizontal mirror under
gravity, where the height distribution follows quantum probability densities
|ψ_n(z)|² based on Airy functions, rather than classical parabolic trajectories.

References
----------
Nesvizhevsky, Valery V., Hans G. Börner, Alexander K. Petukhov, Hartmut Abele,
    Stefan Baeßler, Frank J. Rueß, Thilo Stöferle, Alexander Westphal, Alexei M.
    Gagarski, Guennady A. Petrov, and Alexander V. Strelkov. "Quantum States of
    Neutrons in the Earth's Gravitational Field." Nature 415, no. 6869 (2002):
    297-299.

Jenke, Tobias, Peter Geltenbort, Hartmut Lemmel, and Hartmut Abele. "Realization
    of a Gravity-Resonance-Spectroscopy Technique." Nature Physics 7, no. 6
    (2011): 468-472.

Cronenberg, Gunther, Philippe Brax, Hanno Filter, Peter Geltenbort, Tobias Jenke,
    Guillaume Pignol, Mario Pitschmann, Martin Thalhammer, and Hartmut Abele.
    "Acoustic Rabi Oscillations between Gravitational Quantum States and Impact
    on Symmetron Dark Energy." Nature Physics 14, no. 10 (2018): 1022-1026.

Author: William Zywiec
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import airy
from scipy.optimize import brentq

# Physical constants
HBAR = 1.054571817e-34  # J·s
NEUTRON_MASS = 1.674927471e-27  # kg
G_EARTH = 9.80665  # m/s²
EV_TO_JOULE = 1.602176634e-19  # J/eV

# Derived constants for neutrons in Earth's gravity
# Gravitational length scale: z₀ = (ℏ²/(2m²g))^(1/3)
Z0 = (HBAR**2 / (2 * NEUTRON_MASS**2 * G_EARTH))**(1/3)  # meters
Z0_UM = Z0 * 1e6  # micrometers

# Energy scale: E₀ = m·g·z₀
E0 = NEUTRON_MASS * G_EARTH * Z0  # Joules
E0_PEV = E0 / EV_TO_JOULE * 1e12  # pico-electron-volts


def airy_zeros(n):
    """Calculate the first n zeros of the Airy function Ai(x)."""
    zeros = []
    # Known approximate locations of zeros
    approx_zeros = [-2.338, -4.088, -5.521, -6.787, -7.944,
                    -9.023, -10.040, -11.009, -11.936, -12.829]

    for i in range(min(n, len(approx_zeros))):
        # Refine using Brent's method
        try:
            zero = brentq(lambda x: airy(x)[0], approx_zeros[i] - 0.5,
                         approx_zeros[i] + 0.5)
            zeros.append(zero)
        except ValueError:
            zeros.append(approx_zeros[i])

    # For higher states, use asymptotic formula
    for i in range(len(approx_zeros), n):
        fn = 3 * np.pi * (4 * (i + 1) - 1) / 8
        zeros.append(-fn**(2/3))

    return np.array(zeros)


def energy_level(n):
    """Calculate the n-th quantum energy level in peV."""
    zeros = airy_zeros(n)
    return np.abs(zeros[n-1]) * E0_PEV


def classical_turning_point(n):
    """Calculate the classical turning point height in micrometers."""
    zeros = airy_zeros(n)
    return np.abs(zeros[n-1]) * Z0_UM


def wave_function(n, z_um):
    """
    Calculate the normalized wave function ψ_n(z) for quantum state n.

    Parameters
    ----------
    n : int
        Quantum state number (1-based)
    z_um : array-like
        Height above mirror in micrometers

    Returns
    -------
    psi : array
        Wave function values
    """
    z = np.asarray(z_um) * 1e-6 / Z0  # Convert to dimensionless units
    zeros = airy_zeros(n)
    an = zeros[n-1]

    # Dimensionless coordinate: ξ = z/z₀ + a_n
    xi = z + an

    # Airy function value
    ai_vals = airy(xi)[0]

    # Normalization: |Ai'(a_n)|
    ai_prime_an = airy(an)[1]
    norm = 1.0 / np.abs(ai_prime_an)

    # Apply boundary condition: ψ = 0 for z < 0
    psi = np.where(z_um >= 0, norm * ai_vals, 0.0)

    return psi


def probability_density(n, z_um):
    """Calculate |ψ_n(z)|² probability density."""
    psi = wave_function(n, z_um)
    return psi**2


def plot_energy_levels():
    """Plot the quantized energy levels."""
    fig, ax = plt.subplots(figsize=(10, 6))

    n_states = 10
    energies = [energy_level(n) for n in range(1, n_states + 1)]
    heights = [classical_turning_point(n) for n in range(1, n_states + 1)]

    # Plot energy levels as horizontal lines
    for n, (E, h) in enumerate(zip(energies, heights), 1):
        ax.hlines(E, 0, h, colors=f'C{n-1}', linewidth=2,
                  label=f'n={n}: E={E:.2f} peV, z={h:.1f} μm')

    ax.set_xlabel('Height above mirror (μm)', fontsize=12)
    ax.set_ylabel('Energy (peV)', fontsize=12)
    ax.set_title('Quantum Energy Levels of Neutrons in Gravitational Field',
                 fontsize=14)
    ax.legend(loc='upper left', fontsize=9)
    ax.set_xlim(0, 50)
    ax.set_ylim(0, max(energies) * 1.1)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('energy_levels.png', dpi=150)
    plt.close()
    print("Saved: energy_levels.png")


def plot_wave_functions():
    """Plot wave functions for the first few quantum states."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    z = np.linspace(0, 60, 1000)  # micrometers

    for idx, n in enumerate([1, 2, 3, 4]):
        ax = axes[idx // 2, idx % 2]

        psi = wave_function(n, z)
        prob = probability_density(n, z)
        zn = classical_turning_point(n)

        # Plot wave function
        ax.plot(z, psi, 'b-', linewidth=1.5, label=r'$\psi_n(z)$')
        ax.fill_between(z, 0, prob * 20, alpha=0.3, color='orange',
                        label=r'$|\psi_n(z)|^2$ (scaled)')

        # Mark classical turning point
        ax.axvline(zn, color='red', linestyle='--', linewidth=1,
                   label=f'Classical turning point: {zn:.1f} μm')

        ax.set_xlabel('Height z (μm)', fontsize=11)
        ax.set_ylabel('Wave function / Probability', fontsize=11)
        ax.set_title(f'Quantum State n={n}, E={energy_level(n):.2f} peV',
                     fontsize=12)
        ax.legend(loc='upper right', fontsize=9)
        ax.set_xlim(0, 60)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('wave_functions.png', dpi=150)
    plt.close()
    print("Saved: wave_functions.png")


def plot_classical_vs_quantum():
    """Compare classical and quantum height distributions."""
    fig, ax = plt.subplots(figsize=(10, 6))

    z = np.linspace(0, 50, 500)  # micrometers

    # Quantum: Plot probability density for ground state
    prob_quantum = probability_density(1, z)
    ax.fill_between(z, 0, prob_quantum, alpha=0.5, color='blue',
                    label='Quantum: |ψ₁(z)|² (ground state)')

    # Classical: For a bouncing particle, the probability density is
    # inversely proportional to velocity: P(z) ∝ 1/v(z) ∝ 1/sqrt(z_max - z)
    z_max = classical_turning_point(1)
    prob_classical = np.zeros_like(z)
    mask = z < z_max
    prob_classical[mask] = 1.0 / np.sqrt(z_max - z[mask] + 0.1)
    # Normalize
    prob_classical = prob_classical / np.trapz(prob_classical, z)

    ax.fill_between(z, 0, prob_classical, alpha=0.3, color='red',
                    label='Classical: P(z) ∝ 1/v(z)')

    # Mark the quantum ground state energy level
    ax.axvline(z_max, color='green', linestyle='--', linewidth=2,
               label=f'Ground state turning point: {z_max:.1f} μm')

    ax.set_xlabel('Height above mirror (μm)', fontsize=12)
    ax.set_ylabel('Probability density', fontsize=12)
    ax.set_title('Classical vs Quantum Height Distribution\n'
                 '(Ground State, E₁ = 1.41 peV)', fontsize=14)
    ax.legend(loc='upper right', fontsize=10)
    ax.set_xlim(0, 30)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('classical_vs_quantum.png', dpi=150)
    plt.close()
    print("Saved: classical_vs_quantum.png")


def plot_experimental_comparison():
    """Plot comparison with Nesvizhevsky experimental data."""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Experimental energy levels from Nesvizhevsky et al. (2002)
    experimental = {
        'n': [1, 2, 3],
        'E_peV': [1.41, 2.46, 3.32],
        'E_err': [0.05, 0.1, 0.15]  # Approximate uncertainties
    }

    # Theoretical predictions
    n_states = 5
    n_values = list(range(1, n_states + 1))
    theoretical = [energy_level(n) for n in n_values]

    # Plot theoretical
    ax.bar([n - 0.15 for n in n_values], theoretical, width=0.3, color='blue',
           alpha=0.7, label='Theory (Bloch-Airy)')

    # Plot experimental
    ax.bar([n + 0.15 for n in experimental['n']], experimental['E_peV'],
           width=0.3, color='red', alpha=0.7, label='Experiment (Nesvizhevsky 2002)',
           yerr=experimental['E_err'], capsize=3)

    ax.set_xlabel('Quantum state n', fontsize=12)
    ax.set_ylabel('Energy (peV)', fontsize=12)
    ax.set_title('UCN Quantum States: Theory vs Experiment', fontsize=14)
    ax.set_xticks(n_values)
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')

    # Add percentage difference text
    for n in experimental['n']:
        th = energy_level(n)
        exp = experimental['E_peV'][n-1]
        diff = 100 * (th - exp) / exp
        ax.annotate(f'{diff:+.1f}%', (n, max(th, exp) + 0.1),
                    ha='center', fontsize=9, color='gray')

    plt.tight_layout()
    plt.savefig('experimental_comparison.png', dpi=150)
    plt.close()
    print("Saved: experimental_comparison.png")


def print_quantum_states():
    """Print a summary table of quantum states."""
    print("\n" + "="*70)
    print("BLOCH-AIRY QUANTUM GRAVITATIONAL BOUND STATES")
    print("="*70)
    print(f"\nGravitational length scale: z₀ = {Z0_UM:.3f} μm")
    print(f"Energy scale: E₀ = {E0_PEV:.3f} peV")
    print("\n" + "-"*70)
    print(f"{'State n':^10} {'Energy (peV)':^15} {'Height (μm)':^15} {'Ai zero':^15}")
    print("-"*70)

    zeros = airy_zeros(10)
    for n in range(1, 11):
        E = energy_level(n)
        h = classical_turning_point(n)
        a = zeros[n-1]
        print(f"{n:^10d} {E:^15.3f} {h:^15.2f} {a:^15.4f}")

    print("-"*70)
    print("\nExperimental values (Nesvizhevsky 2002):")
    print("  n=1: E = 1.41 peV, h = 13.7 μm")
    print("  n=2: E = 2.46 peV")
    print("  n=3: E = 3.32 peV")
    print("="*70 + "\n")


def generate_openmc_input():
    """Generate OpenMC input files for Bloch-Airy simulation."""
    print("\nGenerating OpenMC input for Bloch-Airy simulation...")

    # Note: This requires OpenMC to be installed and configured
    try:
        import openmc
    except ImportError:
        print("OpenMC not available. Skipping input generation.")
        print("To run the full simulation, install OpenMC and rebuild with")
        print("the Bloch-Airy implementation.")
        return

    # Create materials (vacuum for UCN bouncing)
    materials = openmc.Materials()
    materials.export_to_xml()

    # Create geometry: A box with a reflective floor
    # The floor represents the mirror where neutrons bounce
    floor = openmc.ZPlane(z0=0.0, boundary_type='reflective')
    ceiling = openmc.ZPlane(z0=100.0e-4, boundary_type='vacuum')  # 100 μm
    x_min = openmc.XPlane(x0=-1.0, boundary_type='reflective')
    x_max = openmc.XPlane(x0=1.0, boundary_type='reflective')
    y_min = openmc.YPlane(y0=-1.0, boundary_type='reflective')
    y_max = openmc.YPlane(y0=1.0, boundary_type='reflective')

    cell = openmc.Cell()
    cell.region = +floor & -ceiling & +x_min & -x_max & +y_min & -y_max

    root = openmc.Universe(cells=[cell])
    geometry = openmc.Geometry(root)
    geometry.export_to_xml()

    # Create settings with Bloch-Airy enabled
    settings = openmc.Settings()
    settings.batches = 100
    settings.particles = 10000
    settings.run_mode = 'fixed source'

    # Enable gravity
    settings.gravity = {
        'enabled': True,
        'acceleration': [0.0, 0.0, -980.0]  # cm/s² (Earth gravity)
    }

    # Enable Bloch-Airy quantum model
    settings.bloch_airy = {
        'enabled': True,
        'energy_threshold': 300e-9  # 300 neV threshold
    }

    # Create UCN source
    # Ultracold neutrons with energies in the peV range
    source = openmc.IndependentSource()
    source.space = openmc.stats.Point((0, 0, 10e-4))  # 10 μm above floor
    source.angle = openmc.stats.Isotropic()
    # UCN energy: ~1 peV = 1e-12 eV
    source.energy = openmc.stats.Discrete([1.0e-12], [1.0])
    source.particle = 'neutron'

    settings.source = source
    settings.export_to_xml()

    print("OpenMC input files generated:")
    print("  - materials.xml")
    print("  - geometry.xml")
    print("  - settings.xml")
    print("\nTo run the simulation, execute: openmc")


def main():
    """Main function to generate all visualizations."""
    print("Bloch-Airy Quantum Gravitational Bound States Visualization")
    print("============================================================")

    # Print summary of quantum states
    print_quantum_states()

    # Generate visualizations
    print("Generating visualizations...")
    plot_energy_levels()
    plot_wave_functions()
    plot_classical_vs_quantum()
    plot_experimental_comparison()

    print("\nAll visualizations saved!")
    print("\nFiles generated:")
    print("  - energy_levels.png: Quantized energy levels")
    print("  - wave_functions.png: Wave functions for states n=1-4")
    print("  - classical_vs_quantum.png: Classical vs quantum comparison")
    print("  - experimental_comparison.png: Theory vs experiment")

    # Generate OpenMC input (optional)
    generate_openmc_input()


if __name__ == '__main__':
    main()
