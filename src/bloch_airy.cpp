#include "openmc/bloch_airy.h"

#include <algorithm>
#include <cmath>

#include "openmc/constants.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"

namespace openmc {

//==============================================================================
// Pre-computed zeros of the Airy function Ai(x)
// These are the first 20 zeros, sufficient for most UCN applications
//==============================================================================

static const double AIRY_ZEROS[] = {
  -2.3381074104597670,  // a_1
  -4.0879494441309710,  // a_2
  -5.5205598280955510,  // a_3
  -6.7867080900819620,  // a_4
  -7.9441335871208230,  // a_5
  -9.0226508533193800,  // a_6
  -10.0401743415580200, // a_7
  -11.0085243037332600, // a_8
  -11.9360155632362000, // a_9
  -12.8287877088894900, // a_10
  -13.6914890352107000, // a_11
  -14.5278299517753800, // a_12
  -15.3407551359779400, // a_13
  -16.1326851569457200, // a_14
  -16.9056339974299000, // a_15
  -17.6613001056970800, // a_16
  -18.4011325992071500, // a_17
  -19.1263804742469400, // a_18
  -19.8381298917214600, // a_19
  -20.5373329076775400  // a_20
};

static const int NUM_PRECOMPUTED_ZEROS = 20;

//==============================================================================
// Airy function implementation using polynomial/asymptotic approximations
//==============================================================================

double airy_ai(double x)
{
  // For x <= 0: Use power series expansion
  // For x > 0: Use asymptotic expansion

  if (x < -5.0) {
    // Asymptotic expansion for large negative x
    double z = std::pow(-x, 1.5) * 2.0 / 3.0;
    double sqrtx = std::sqrt(-x);
    double c = std::cos(z - M_PI / 4.0);
    double s = std::sin(z - M_PI / 4.0);

    // Leading term of asymptotic expansion
    double ai = (c + s / (48.0 * z)) / (std::sqrt(M_PI) * std::pow(-x, 0.25));
    return ai;
  } else if (x < 0.0) {
    // Power series for small negative x
    // Ai(x) = c1 * f(x) - c2 * g(x)
    // where c1 = Ai(0) ≈ 0.355028, c2 = -Ai'(0) ≈ 0.258819
    double c1 = 0.35502805388781723;
    double c2 = 0.25881940379280680;

    double x3 = x * x * x;
    double f = 1.0;
    double g = x;
    double term_f = 1.0;
    double term_g = x;

    for (int k = 1; k <= 30; ++k) {
      term_f *= x3 / ((3 * k - 1) * (3 * k));
      term_g *= x3 / ((3 * k) * (3 * k + 1));
      f += term_f;
      g += term_g;
      if (std::abs(term_f) < 1e-15 && std::abs(term_g) < 1e-15) break;
    }

    return c1 * f - c2 * g;
  } else if (x < 2.0) {
    // Power series for small positive x
    double c1 = 0.35502805388781723;
    double c2 = 0.25881940379280680;

    double x3 = x * x * x;
    double f = 1.0;
    double g = x;
    double term_f = 1.0;
    double term_g = x;

    for (int k = 1; k <= 30; ++k) {
      term_f *= x3 / ((3 * k - 1) * (3 * k));
      term_g *= x3 / ((3 * k) * (3 * k + 1));
      f += term_f;
      g += term_g;
      if (std::abs(term_f) < 1e-15 && std::abs(term_g) < 1e-15) break;
    }

    return c1 * f - c2 * g;
  } else {
    // Asymptotic expansion for large positive x
    double z = std::pow(x, 1.5) * 2.0 / 3.0;
    double ai = std::exp(-z) / (2.0 * std::sqrt(M_PI) * std::pow(x, 0.25));

    // Add correction terms
    double u = 1.0 / z;
    double sum = 1.0 - 5.0 * u / 72.0;

    return ai * sum;
  }
}

double airy_ai_prime(double x)
{
  // Numerical derivative using central difference
  double h = 1e-6;
  return (airy_ai(x + h) - airy_ai(x - h)) / (2.0 * h);
}

double airy_zero(int n)
{
  if (n < 1) {
    return AIRY_ZEROS[0];
  }

  if (n <= NUM_PRECOMPUTED_ZEROS) {
    return AIRY_ZEROS[n - 1];
  }

  // Asymptotic approximation for large n:
  // a_n ≈ -f_n where f_n = (3π(4n-1)/8)^(2/3)
  double fn = 3.0 * M_PI * (4 * n - 1) / 8.0;
  return -std::pow(fn, 2.0 / 3.0);
}

//==============================================================================
// Gravitational quantum mechanics functions
//==============================================================================

double gravitational_length_scale(double mass_ev, double g_cm_s2)
{
  // z₀ = (ℏ²c²/(2m²g))^(1/3)
  // Note: We need to be careful with units here
  // ℏ in eV·s, m in eV/c², g in cm/s², result in cm

  double hbar_c = HBAR_EV_S * C_CM_S;  // eV·cm
  double numerator = hbar_c * hbar_c;
  double denominator = 2.0 * mass_ev * mass_ev * g_cm_s2;

  return std::pow(numerator / denominator, 1.0 / 3.0);
}

double gravitational_energy_scale(double mass_ev, double g_cm_s2)
{
  // E₀ = m·g·z₀ / c²
  double z0 = gravitational_length_scale(mass_ev, g_cm_s2);
  return mass_ev * g_cm_s2 * z0 / (C_CM_S * C_CM_S);
}

double quantum_energy_level(int n, double mass_ev, double g_cm_s2)
{
  double E0 = gravitational_energy_scale(mass_ev, g_cm_s2);
  double an = airy_zero(n);
  return std::abs(an) * E0;
}

double classical_turning_point(int n, double mass_ev, double g_cm_s2)
{
  double z0 = gravitational_length_scale(mass_ev, g_cm_s2);
  double an = airy_zero(n);
  return std::abs(an) * z0;
}

double wave_function(int n, double z, double z0)
{
  if (z < 0.0) {
    return 0.0;  // Wave function is zero below the mirror
  }

  double an = airy_zero(n);
  double zn = std::abs(an) * z0;  // Classical turning point

  // Dimensionless coordinate: ξ = (z - z_n) / z₀ = z/z₀ + a_n
  double xi = z / z0 + an;

  // Normalization constant: C_n = 1 / (z₀^(1/2) * |Ai'(a_n)|)
  double ai_prime_an = airy_ai_prime(an);
  double Cn = 1.0 / (std::sqrt(z0) * std::abs(ai_prime_an));

  return Cn * airy_ai(xi);
}

double probability_density(int n, double z, double z0)
{
  double psi = wave_function(n, z, z0);
  return psi * psi;
}

double sample_quantum_height(int n, double z0, uint64_t* seed)
{
  // Use rejection sampling from the probability distribution |ψ_n(z)|²
  // The distribution is bounded by the classical turning point

  double zn = classical_turning_point(n, NEUTRON_MASS_EV,
    std::sqrt(settings::gravity_accel[0] * settings::gravity_accel[0] +
              settings::gravity_accel[1] * settings::gravity_accel[1] +
              settings::gravity_accel[2] * settings::gravity_accel[2]));

  // Maximum height to sample (extend slightly beyond classical turning point)
  double z_max = zn * 1.5;

  // Find maximum of probability density for rejection sampling
  // The maximum is near z = 0 for all states
  double p_max = probability_density(n, 0.0, z0) * 1.1;  // Add safety margin

  // Rejection sampling
  int max_iterations = 10000;
  for (int i = 0; i < max_iterations; ++i) {
    double z_try = z_max * prn(seed);
    double p_try = probability_density(n, z_try, z0);
    double accept = p_max * prn(seed);

    if (accept < p_try) {
      return z_try;
    }
  }

  // Fallback: return classical turning point
  return zn * 0.5;
}

int determine_quantum_state(double energy, double mass_ev, double g_cm_s2)
{
  // Find quantum state n such that E_n is closest to the given energy
  double E0 = gravitational_energy_scale(mass_ev, g_cm_s2);

  // E_n = |a_n| * E0, so |a_n| = E / E0
  double target_an = energy / E0;

  // Search through quantum states
  int n = 1;
  while (n < NUM_PRECOMPUTED_ZEROS) {
    double an_curr = std::abs(airy_zero(n));
    double an_next = std::abs(airy_zero(n + 1));

    if (target_an < (an_curr + an_next) / 2.0) {
      break;
    }
    ++n;
  }

  return std::max(1, n);
}

bool should_apply_bloch_airy(int particle_type, double energy)
{
  // Check if Bloch-Airy is enabled
  if (!settings::bloch_airy_enabled) {
    return false;
  }

  // Check if gravity is enabled (required for Bloch-Airy)
  if (!settings::gravity_enabled) {
    return false;
  }

  // Check if particle is a neutron (particle_type == 0 for neutrons in OpenMC)
  if (particle_type != 0) {
    return false;
  }

  // Check if energy is below threshold
  if (energy > settings::bloch_airy_energy_threshold) {
    return false;
  }

  return true;
}

} // namespace openmc
