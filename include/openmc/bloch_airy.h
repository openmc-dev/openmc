//! \file bloch_airy.h
//! Bloch-Airy quantum gravitational bound states for ultracold neutrons.
//!
//! This module provides functions for calculating quantum mechanical effects
//! for ultracold neutrons (UCN) bouncing in a gravitational potential. The
//! solutions to the Schrödinger equation in a linear potential are Airy
//! functions, with quantized energy levels determined by the zeros of Ai(x).

#ifndef OPENMC_BLOCH_AIRY_H
#define OPENMC_BLOCH_AIRY_H

#include <cstdint>

namespace openmc {

//==============================================================================
// Physical constants for Bloch-Airy calculations
//==============================================================================

constexpr double HBAR_EV_S = 6.582119569e-16;    //!< Reduced Planck constant [eV·s]
constexpr double NEUTRON_MASS_EV = 939.565346e6; //!< Neutron mass [eV/c²]
constexpr double C_CM_S = 2.99792458e10;         //!< Speed of light [cm/s]

//==============================================================================
//! Calculate the Airy function Ai(x).
//!
//! Uses a polynomial/asymptotic approximation for the Airy function of the
//! first kind, which appears in the wave functions of particles in a linear
//! potential (gravity).
//!
//! \param x The argument of the Airy function
//! \return The value of Ai(x)
//==============================================================================

double airy_ai(double x);

//==============================================================================
//! Calculate the derivative of the Airy function Ai'(x).
//!
//! \param x The argument of the Airy function
//! \return The value of Ai'(x)
//==============================================================================

double airy_ai_prime(double x);

//==============================================================================
//! Get the n-th zero of the Airy function Ai(x).
//!
//! The zeros of Ai(x) determine the quantized energy levels of particles in
//! a linear gravitational potential. The first few zeros are:
//!   a_1 ≈ -2.338, a_2 ≈ -4.088, a_3 ≈ -5.521, ...
//!
//! \param n The index of the zero (1-based: n=1 is the first zero)
//! \return The n-th zero of Ai(x) (negative value)
//==============================================================================

double airy_zero(int n);

//==============================================================================
//! Calculate the gravitational length scale z₀.
//!
//! z₀ = (ℏ²/(2m²g))^(1/3)
//!
//! For neutrons in Earth's gravity, z₀ ≈ 5.87 μm.
//!
//! \param mass_ev Particle mass in eV/c²
//! \param g_cm_s2 Gravitational acceleration magnitude in cm/s²
//! \return The gravitational length scale z₀ in cm
//==============================================================================

double gravitational_length_scale(double mass_ev, double g_cm_s2);

//==============================================================================
//! Calculate the gravitational energy scale E₀.
//!
//! E₀ = m·g·z₀
//!
//! For neutrons in Earth's gravity, E₀ ≈ 0.602 peV.
//!
//! \param mass_ev Particle mass in eV/c²
//! \param g_cm_s2 Gravitational acceleration magnitude in cm/s²
//! \return The gravitational energy scale E₀ in eV
//==============================================================================

double gravitational_energy_scale(double mass_ev, double g_cm_s2);

//==============================================================================
//! Calculate the n-th quantized energy level.
//!
//! E_n = |a_n| × E₀
//!
//! where a_n is the n-th zero of Ai(x).
//!
//! \param n The quantum number (1-based)
//! \param mass_ev Particle mass in eV/c²
//! \param g_cm_s2 Gravitational acceleration magnitude in cm/s²
//! \return The n-th energy eigenvalue in eV
//==============================================================================

double quantum_energy_level(int n, double mass_ev, double g_cm_s2);

//==============================================================================
//! Calculate the classical turning point height for state n.
//!
//! z_n = |a_n| × z₀
//!
//! This is the maximum height a classical particle with energy E_n would reach.
//!
//! \param n The quantum number (1-based)
//! \param mass_ev Particle mass in eV/c²
//! \param g_cm_s2 Gravitational acceleration magnitude in cm/s²
//! \return The classical turning point height in cm
//==============================================================================

double classical_turning_point(int n, double mass_ev, double g_cm_s2);

//==============================================================================
//! Calculate the normalized wave function ψ_n(z) for quantum state n.
//!
//! ψ_n(z) = C_n × Ai((z - z_n)/z₀)
//!
//! where C_n is the normalization constant.
//!
//! \param n The quantum number (1-based)
//! \param z Height above the mirror in cm
//! \param z0 Gravitational length scale in cm
//! \return The wave function value (normalized)
//==============================================================================

double wave_function(int n, double z, double z0);

//==============================================================================
//! Calculate the probability density |ψ_n(z)|² for quantum state n.
//!
//! \param n The quantum number (1-based)
//! \param z Height above the mirror in cm
//! \param z0 Gravitational length scale in cm
//! \return The probability density
//==============================================================================

double probability_density(int n, double z, double z0);

//==============================================================================
//! Sample a height from the quantum probability distribution |ψ_n(z)|².
//!
//! Uses rejection sampling to draw a random height from the probability
//! distribution of quantum state n.
//!
//! \param n The quantum number (1-based)
//! \param z0 Gravitational length scale in cm
//! \param seed Pointer to the random number generator seed
//! \return A randomly sampled height in cm
//==============================================================================

double sample_quantum_height(int n, double z0, uint64_t* seed);

//==============================================================================
//! Determine the quantum state number for a given neutron energy.
//!
//! Finds the quantum state n whose energy E_n is closest to the given energy.
//!
//! \param energy Neutron kinetic energy in eV
//! \param mass_ev Particle mass in eV/c²
//! \param g_cm_s2 Gravitational acceleration magnitude in cm/s²
//! \return The quantum state number (1-based)
//==============================================================================

int determine_quantum_state(double energy, double mass_ev, double g_cm_s2);

//==============================================================================
//! Check if Bloch-Airy quantum effects should be applied.
//!
//! Returns true if:
//!   1. Bloch-Airy is enabled in settings
//!   2. The particle is a neutron
//!   3. The neutron energy is below the threshold
//!
//! \param particle_type The type of particle
//! \param energy The particle energy in eV
//! \return True if quantum effects should be applied
//==============================================================================

bool should_apply_bloch_airy(int particle_type, double energy);

} // namespace openmc

#endif // OPENMC_BLOCH_AIRY_H
