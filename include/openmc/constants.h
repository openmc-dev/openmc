//! \file constants.h
//! A collection of constants

#ifndef OPENMC_CONSTANTS_H
#define OPENMC_CONSTANTS_H

#include <cmath>
#include <cstdint>
#include <limits>

#include "openmc/array.h"
#include "openmc/vector.h"
#include "openmc/version.h"

namespace openmc {

using double_2dvec = vector<vector<double>>;
using double_3dvec = vector<vector<vector<double>>>;
using double_4dvec = vector<vector<vector<vector<double>>>>;

// ============================================================================
// VERSIONING NUMBERS

// HDF5 data format
constexpr int HDF5_VERSION[] {3, 0};

// Version numbers for binary files
constexpr array<int, 2> VERSION_STATEPOINT {18, 1};
constexpr array<int, 2> VERSION_PARTICLE_RESTART {2, 0};
constexpr array<int, 2> VERSION_TRACK {3, 0};
constexpr array<int, 2> VERSION_SUMMARY {6, 0};
constexpr array<int, 2> VERSION_VOLUME {1, 0};
constexpr array<int, 2> VERSION_VOXEL {2, 0};
constexpr array<int, 2> VERSION_MGXS_LIBRARY {1, 0};
constexpr array<int, 2> VERSION_PROPERTIES {1, 0};
constexpr array<int, 2> VERSION_WEIGHT_WINDOWS {1, 0};

// ============================================================================
// ADJUSTABLE PARAMETERS

// NOTE: This is the only section of the constants module that should ever be
// adjusted. Modifying constants in other sections may cause the code to fail.

// Significance level for confidence intervals
constexpr double CONFIDENCE_LEVEL {0.95};

// Used for surface current tallies
constexpr double TINY_BIT {1e-8};

// User for precision in geometry
constexpr double FP_PRECISION {1e-14};
constexpr double FP_REL_PRECISION {1e-5};
constexpr double FP_COINCIDENT {1e-12};

// Coincidence tolerances
constexpr double TORUS_TOL {1e-10};
constexpr double RADIAL_MESH_TOL {1e-10};

// Maximum number of random samples per history
constexpr int MAX_SAMPLE {100000};

// ============================================================================
// MATH AND PHYSICAL CONSTANTS

// TODO: replace with <numbers> when we go for C++20
constexpr double PI {3.141592653589793238462643383279502884L};
const double SQRT_PI {std::sqrt(PI)};
constexpr double INFTY {std::numeric_limits<double>::max()};

// Values here are from the Committee on Data for Science and Technology
// (CODATA) 2018 recommendation (https://physics.nist.gov/cuu/Constants/).

// Physical constants
constexpr double MASS_NEUTRON {1.00866491595}; // mass of a neutron in amu
constexpr double MASS_NEUTRON_EV {
  939.56542052e6};                             // mass of a neutron in eV/c^2
constexpr double MASS_PROTON {1.007276466621}; // mass of a proton in amu
constexpr double MASS_ELECTRON_EV {
  0.51099895000e6}; // electron mass energy equivalent in eV/c^2
constexpr double FINE_STRUCTURE {
  137.035999084}; // inverse fine structure constant
constexpr double PLANCK_C {
  1.2398419839593942e4}; // Planck's constant times c in eV-Angstroms
constexpr double AMU {1.66053906660e-27};      // 1 amu in kg
constexpr double C_LIGHT {2.99792458e10};      // speed of light in cm/s
constexpr double N_AVOGADRO {0.602214076};     // Avogadro's number in 10^24/mol
constexpr double K_BOLTZMANN {8.617333262e-5}; // Boltzmann constant in eV/K

// Electron subshell labels
constexpr array<const char*, 39> SUBSHELLS = {"K", "L1", "L2", "L3", "M1", "M2",
  "M3", "M4", "M5", "N1", "N2", "N3", "N4", "N5", "N6", "N7", "O1", "O2", "O3",
  "O4", "O5", "O6", "O7", "O8", "O9", "P1", "P2", "P3", "P4", "P5", "P6", "P7",
  "P8", "P9", "P10", "P11", "Q1", "Q2", "Q3"};

// Void material and nuclide
// TODO: refactor and remove
constexpr int MATERIAL_VOID {-1};
constexpr int NUCLIDE_NONE {-1};

// ============================================================================
// CROSS SECTION RELATED CONSTANTS

// Temperature treatment method
enum class TemperatureMethod { NEAREST, INTERPOLATION };

// Reaction types
enum ReactionType {
  REACTION_NONE = 0,
  TOTAL_XS = 1,
  ELASTIC = 2,
  N_NONELASTIC = 3,
  N_LEVEL = 4,
  MISC = 5,
  N_2ND = 11,
  N_2N = 16,
  N_3N = 17,
  N_FISSION = 18,
  N_F = 19,
  N_NF = 20,
  N_2NF = 21,
  N_NA = 22,
  N_N3A = 23,
  N_2NA = 24,
  N_3NA = 25,
  N_NP = 28,
  N_N2A = 29,
  N_2N2A = 30,
  N_ND = 32,
  N_NT = 33,
  N_N3HE = 34,
  N_ND2A = 35,
  N_NT2A = 36,
  N_4N = 37,
  N_3NF = 38,
  N_2NP = 41,
  N_3NP = 42,
  N_N2P = 44,
  N_NPA = 45,
  N_N1 = 51,
  N_N40 = 90,
  N_NC = 91,
  N_DISAPPEAR = 101,
  N_GAMMA = 102,
  N_P = 103,
  N_D = 104,
  N_T = 105,
  N_3HE = 106,
  N_A = 107,
  N_2A = 108,
  N_3A = 109,
  N_2P = 111,
  N_PA = 112,
  N_T2A = 113,
  N_D2A = 114,
  N_PD = 115,
  N_PT = 116,
  N_DA = 117,
  N_5N = 152,
  N_6N = 153,
  N_2NT = 154,
  N_TA = 155,
  N_4NP = 156,
  N_3ND = 157,
  N_NDA = 158,
  N_2NPA = 159,
  N_7N = 160,
  N_8N = 161,
  N_5NP = 162,
  N_6NP = 163,
  N_7NP = 164,
  N_4NA = 165,
  N_5NA = 166,
  N_6NA = 167,
  N_7NA = 168,
  N_4ND = 169,
  N_5ND = 170,
  N_6ND = 171,
  N_3NT = 172,
  N_4NT = 173,
  N_5NT = 174,
  N_6NT = 175,
  N_2N3HE = 176,
  N_3N3HE = 177,
  N_4N3HE = 178,
  N_3N2P = 179,
  N_3N2A = 180,
  N_3NPA = 181,
  N_DT = 182,
  N_NPD = 183,
  N_NPT = 184,
  N_NDT = 185,
  N_NP3HE = 186,
  N_ND3HE = 187,
  N_NT3HE = 188,
  N_NTA = 189,
  N_2N2P = 190,
  N_P3HE = 191,
  N_D3HE = 192,
  N_3HEA = 193,
  N_4N2P = 194,
  N_4N2A = 195,
  N_4NPA = 196,
  N_3P = 197,
  N_N3P = 198,
  N_3N2PA = 199,
  N_5N2P = 200,
  N_XP = 203,
  N_XD = 204,
  N_XT = 205,
  N_X3HE = 206,
  N_XA = 207,
  HEATING = 301,
  DAMAGE_ENERGY = 444,
  COHERENT = 502,
  INCOHERENT = 504,
  PAIR_PROD_ELEC = 515,
  PAIR_PROD = 516,
  PAIR_PROD_NUC = 517,
  PHOTOELECTRIC = 522,
  N_P0 = 600,
  N_PC = 649,
  N_D0 = 650,
  N_DC = 699,
  N_T0 = 700,
  N_TC = 749,
  N_3HE0 = 750,
  N_3HEC = 799,
  N_A0 = 800,
  N_AC = 849,
  N_2N0 = 875,
  N_2NC = 891,
  HEATING_LOCAL = 901
};

constexpr array<int, 6> DEPLETION_RX {N_GAMMA, N_P, N_A, N_2N, N_3N, N_4N};

// Resonance elastic scattering methods
enum class ResScatMethod {
  rvs,  // Relative velocity sampling
  dbrc, // Doppler broadening rejection correction
  cxs   // Constant cross section
};

enum class ElectronTreatment {
  LED, // Local Energy Deposition
  TTB  // Thick Target Bremsstrahlung
};

// ============================================================================
// MULTIGROUP RELATED

// Flag to denote this was a macroscopic data object
constexpr double MACROSCOPIC_AWR {-2.};

// Number of mu bins to use when converting Legendres to tabular type
constexpr int DEFAULT_NMU {33};

// Mgxs::get_xs enumerated types
enum class MgxsType {
  TOTAL,
  ABSORPTION,
  INVERSE_VELOCITY,
  DECAY_RATE,
  NU_SCATTER,
  SCATTER,
  NU_SCATTER_FMU,
  SCATTER_FMU,
  FISSION,
  KAPPA_FISSION,
  PROMPT_NU_FISSION,
  DELAYED_NU_FISSION,
  NU_FISSION,
  CHI_PROMPT,
  CHI_DELAYED
};

// ============================================================================
// TALLY-RELATED CONSTANTS

enum class TallyResult { VALUE, SUM, SUM_SQ, SIZE };

enum class TallyType { VOLUME, MESH_SURFACE, SURFACE, PULSE_HEIGHT };

enum class TallyEstimator { ANALOG, TRACKLENGTH, COLLISION };

enum class TallyEvent { SURFACE, LATTICE, KILL, SCATTER, ABSORB };

// Tally score type -- if you change these, make sure you also update the
// _SCORES dictionary in openmc/capi/tally.py
//
// These are kept as a normal enum and made negative, since variables which
// store one of these enum values usually also may be responsible for storing
// MT numbers from the long enum above.
enum TallyScore {
  SCORE_FLUX = -1,                // flux
  SCORE_TOTAL = -2,               // total reaction rate
  SCORE_SCATTER = -3,             // scattering rate
  SCORE_NU_SCATTER = -4,          // scattering production rate
  SCORE_ABSORPTION = -5,          // absorption rate
  SCORE_FISSION = -6,             // fission rate
  SCORE_NU_FISSION = -7,          // neutron production rate
  SCORE_KAPPA_FISSION = -8,       // fission energy production rate
  SCORE_CURRENT = -9,             // current
  SCORE_EVENTS = -10,             // number of events
  SCORE_DELAYED_NU_FISSION = -11, // delayed neutron production rate
  SCORE_PROMPT_NU_FISSION = -12,  // prompt neutron production rate
  SCORE_INVERSE_VELOCITY = -13,   // flux-weighted inverse velocity
  SCORE_FISS_Q_PROMPT = -14,      // prompt fission Q-value
  SCORE_FISS_Q_RECOV = -15,       // recoverable fission Q-value
  SCORE_DECAY_RATE = -16,         // delayed neutron precursor decay rate
  SCORE_PULSE_HEIGHT = -17        // pulse-height
};

// Global tally parameters
constexpr int N_GLOBAL_TALLIES {4};
enum class GlobalTally { K_COLLISION, K_ABSORPTION, K_TRACKLENGTH, LEAKAGE };

// Miscellaneous
constexpr int C_NONE {-1};

// Interpolation rules
enum class Interpolation {
  histogram = 1,
  lin_lin = 2,
  lin_log = 3,
  log_lin = 4,
  log_log = 5,
  // skip 6 b/c ENDF-6 reserves this value for
  // "special one-dimensional interpolation law"
  quadratic = 7,
  cubic = 8
};

enum class RunMode {
  UNSET, // default value, OpenMC throws error if left to this
  FIXED_SOURCE,
  EIGENVALUE,
  PLOTTING,
  PARTICLE,
  VOLUME
};

enum class SolverType { MONTE_CARLO, RANDOM_RAY };

enum class RandomRaySourceShape { FLAT, LINEAR, LINEAR_XY };

//==============================================================================
// Geometry Constants

enum class GeometryType { CSG, DAG };

//==============================================================================
// Shared Libraries on Windows need the functions/classes/variables to be
// exposed to be declared with this special keyword "__declspec( dllexport )".
// We use a macro for this, so that it is only exported on Windows.

#if _MSC_VER && OPENMC_WIN_COMPILE
#define DllExport __declspec(dllexport)
#elif _MSC_VER
#define DllExport __declspec(dllimport)
#else
#define DllExport
#endif

} // namespace openmc

#endif // OPENMC_CONSTANTS_H
