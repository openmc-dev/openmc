//! \file constants.h
//! A collection of constants

#ifndef OPENMC_CONSTANTS_H
#define OPENMC_CONSTANTS_H

#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace openmc {

using double_2dvec = std::vector<std::vector<double>>;
using double_3dvec = std::vector<std::vector<std::vector<double>>>;
using double_4dvec = std::vector<std::vector<std::vector<std::vector<double>>>>;

// ============================================================================
// VERSIONING NUMBERS

// OpenMC major, minor, and release numbers
constexpr int VERSION_MAJOR {0};
constexpr int VERSION_MINOR {10};
constexpr int VERSION_RELEASE {0};
constexpr std::array<int, 3> VERSION {VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE};

// HDF5 data format
constexpr int HDF5_VERSION[] {2, 0};

// Version numbers for binary files
constexpr std::array<int, 2> VERSION_STATEPOINT {17, 0};
constexpr std::array<int, 2> VERSION_PARTICLE_RESTART {2, 0};
constexpr std::array<int, 2> VERSION_TRACK {2, 0};
constexpr std::array<int, 2> VERSION_SUMMARY {6, 0};
constexpr std::array<int, 2> VERSION_VOLUME {1, 0};
constexpr std::array<int, 2> VERSION_VOXEL {2, 0};
constexpr std::array<int, 2> VERSION_MGXS_LIBRARY {1, 0};
constexpr char VERSION_MULTIPOLE[] {"v0.2"};

// ============================================================================
// ADJUSTABLE PARAMETERS

// NOTE: This is the only section of the constants module that should ever be
// adjusted. Modifying constants in other sections may cause the code to fail.

// Monoatomic ideal-gas scattering treatment threshold
constexpr double FREE_GAS_THRESHOLD {400.0};

// Significance level for confidence intervals
constexpr double CONFIDENCE_LEVEL {0.95};

// Used for surface current tallies
constexpr double TINY_BIT {1e-8};

// User for precision in geometry
constexpr double FP_PRECISION {1e-14};
constexpr double FP_REL_PRECISION {1e-5};
constexpr double FP_COINCIDENT {1e-12};

// Maximum number of collisions/crossings
constexpr int MAX_EVENTS {1000000};
constexpr int MAX_SAMPLE {100000};

// Maximum number of words in a single line, length of line, and length of
// single word
constexpr int MAX_LINE_LEN {250};
constexpr int MAX_WORD_LEN {150};

// Maximum number of external source spatial resamples to encounter before an
// error is thrown.
constexpr int EXTSRC_REJECT_THRESHOLD {10000};
constexpr double EXTSRC_REJECT_FRACTION {0.05};

// ============================================================================
// MATH AND PHYSICAL CONSTANTS

// Values here are from the Committee on Data for Science and Technology
// (CODATA) 2014 recommendation (doi:10.1103/RevModPhys.88.035009).

// TODO: cmath::M_PI has 3 more digits precision than the Fortran constant we
// use so for now we will reuse the Fortran constant until we are OK with
// modifying test results
constexpr double PI {3.1415926535898};
const double SQRT_PI {std::sqrt(PI)};
constexpr double INFTY {std::numeric_limits<double>::max()};

// Physical constants
constexpr double MASS_NEUTRON     {1.00866491588}; // mass of a neutron in amu
constexpr double MASS_NEUTRON_EV  {939.5654133e6}; // mass of a neutron in eV/c^2
constexpr double MASS_PROTON      {1.007276466879}; // mass of a proton in amu
constexpr double MASS_ELECTRON_EV {0.5109989461e6}; // electron mass energy equivalent in eV/c^2
constexpr double FINE_STRUCTURE   {137.035999139}; // inverse fine structure constant
constexpr double PLANCK_C         {1.2398419739062977e4}; // Planck's constant times c in eV-Angstroms
constexpr double AMU              {1.660539040e-27}; // 1 amu in kg
constexpr double C_LIGHT          {2.99792458e8}; // speed of light in m/s
constexpr double N_AVOGADRO       {0.6022140857}; // Avogadro's number in 10^24/mol
constexpr double K_BOLTZMANN      {8.6173303e-5}; // Boltzmann constant in eV/K

// Electron subshell labels
constexpr std::array<const char*, 39> SUBSHELLS =  {
  "K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5",
  "N1", "N2", "N3", "N4", "N5", "N6", "N7", "O1", "O2",
  "O3", "O4", "O5", "O6", "O7", "O8", "O9", "P1", "P2",
  "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11",
  "Q1", "Q2", "Q3"
};

// Void material and nuclide
// TODO: refactor and remove
constexpr int MATERIAL_VOID {-1};
constexpr int NUCLIDE_NONE  {-1};

// ============================================================================
// CROSS SECTION RELATED CONSTANTS

// Angular distribution type
// TODO: Convert to enum
constexpr int ANGLE_ISOTROPIC {1};
constexpr int ANGLE_32_EQUI {2};
constexpr int ANGLE_TABULAR {3};
constexpr int ANGLE_LEGENDRE {4};
constexpr int ANGLE_HISTOGRAM {5};

// Temperature treatment method
// TODO: Convert to enum?
constexpr int TEMPERATURE_NEAREST {1};
constexpr int TEMPERATURE_INTERPOLATION {2};

// Reaction types
// TODO: Convert to enum
constexpr int REACTION_NONE {0};
constexpr int TOTAL_XS {1};
constexpr int ELASTIC {2};
constexpr int N_NONELASTIC {3};
constexpr int N_LEVEL {4};
constexpr int MISC    {5};
constexpr int N_2ND   {11};
constexpr int N_2N    {16};
constexpr int N_3N    {17};
constexpr int N_FISSION {18};
constexpr int N_F     {19};
constexpr int N_NF    {20};
constexpr int N_2NF   {21};
constexpr int N_NA    {22};
constexpr int N_N3A   {23};
constexpr int N_2NA   {24};
constexpr int N_3NA   {25};
constexpr int N_NP    {28};
constexpr int N_N2A   {29};
constexpr int N_2N2A  {30};
constexpr int N_ND    {32};
constexpr int N_NT    {33};
constexpr int N_N3HE  {34};
constexpr int N_ND2A  {35};
constexpr int N_NT2A  {36};
constexpr int N_4N    {37};
constexpr int N_3NF   {38};
constexpr int N_2NP   {41};
constexpr int N_3NP   {42};
constexpr int N_N2P   {44};
constexpr int N_NPA   {45};
constexpr int N_N1    {51};
constexpr int N_N40   {90};
constexpr int N_NC    {91};
constexpr int N_DISAPPEAR {101};
constexpr int N_GAMMA {102};
constexpr int N_P     {103};
constexpr int N_D     {104};
constexpr int N_T     {105};
constexpr int N_3HE   {106};
constexpr int N_A     {107};
constexpr int N_2A    {108};
constexpr int N_3A    {109};
constexpr int N_2P    {111};
constexpr int N_PA    {112};
constexpr int N_T2A   {113};
constexpr int N_D2A   {114};
constexpr int N_PD    {115};
constexpr int N_PT    {116};
constexpr int N_DA    {117};
constexpr int N_5N    {152};
constexpr int N_6N    {153};
constexpr int N_2NT   {154};
constexpr int N_TA    {155};
constexpr int N_4NP   {156};
constexpr int N_3ND   {157};
constexpr int N_NDA   {158};
constexpr int N_2NPA  {159};
constexpr int N_7N    {160};
constexpr int N_8N    {161};
constexpr int N_5NP   {162};
constexpr int N_6NP   {163};
constexpr int N_7NP   {164};
constexpr int N_4NA   {165};
constexpr int N_5NA   {166};
constexpr int N_6NA   {167};
constexpr int N_7NA   {168};
constexpr int N_4ND   {169};
constexpr int N_5ND   {170};
constexpr int N_6ND   {171};
constexpr int N_3NT   {172};
constexpr int N_4NT   {173};
constexpr int N_5NT   {174};
constexpr int N_6NT   {175};
constexpr int N_2N3HE {176};
constexpr int N_3N3HE {177};
constexpr int N_4N3HE {178};
constexpr int N_3N2P  {179};
constexpr int N_3N3A  {180};
constexpr int N_3NPA  {181};
constexpr int N_DT    {182};
constexpr int N_NPD   {183};
constexpr int N_NPT   {184};
constexpr int N_NDT   {185};
constexpr int N_NP3HE {186};
constexpr int N_ND3HE {187};
constexpr int N_NT3HE {188};
constexpr int N_NTA   {189};
constexpr int N_2N2P  {190};
constexpr int N_P3HE  {191};
constexpr int N_D3HE  {192};
constexpr int N_3HEA  {193};
constexpr int N_4N2P  {194};
constexpr int N_4N2A  {195};
constexpr int N_4NPA  {196};
constexpr int N_3P    {197};
constexpr int N_N3P   {198};
constexpr int N_3N2PA {199};
constexpr int N_5N2P  {200};
constexpr int N_XP    {203};
constexpr int N_XD    {204};
constexpr int N_XT    {205};
constexpr int N_X3HE  {206};
constexpr int N_XA    {207};
constexpr int NEUTRON_HEATING {301};
constexpr int DAMAGE_ENERGY {444};
constexpr int COHERENT {502};
constexpr int INCOHERENT {504};
constexpr int PAIR_PROD_ELEC {515};
constexpr int PAIR_PROD {516};
constexpr int PAIR_PROD_NUC {517};
constexpr int PHOTOELECTRIC {522};
constexpr int N_P0    {600};
constexpr int N_PC    {649};
constexpr int N_D0    {650};
constexpr int N_DC    {699};
constexpr int N_T0    {700};
constexpr int N_TC    {749};
constexpr int N_3HE0  {750};
constexpr int N_3HEC  {799};
constexpr int N_A0    {800};
constexpr int N_AC    {849};
constexpr int N_2N0   {875};
constexpr int N_2NC   {891};

constexpr std::array<int, 6> DEPLETION_RX {N_GAMMA, N_P, N_A, N_2N, N_3N, N_4N};

// Fission neutron emission (nu) type
constexpr int NU_NONE       {0}; // No nu values (non-fissionable)
constexpr int NU_POLYNOMIAL {1}; // Nu values given by polynomial
constexpr int NU_TABULAR    {2}; // Nu values given by tabular distribution

// Library types
constexpr int LIBRARY_NEUTRON {1};
constexpr int LIBRARY_THERMAL {2};
constexpr int LIBRARY_PHOTON {3};
constexpr int LIBRARY_MULTIGROUP {4};

// Probability table parameters
constexpr int URR_CUM_PROB {0};
constexpr int URR_TOTAL    {1};
constexpr int URR_ELASTIC  {2};
constexpr int URR_FISSION  {3};
constexpr int URR_N_GAMMA  {4};
constexpr int URR_HEATING  {5};

// Maximum number of partial fission reactions
constexpr int PARTIAL_FISSION_MAX {4};

// Resonance elastic scattering methods
// TODO: Convert to enum
enum class ResScatMethod {
  rvs, // Relative velocity sampling
  dbrc, // Doppler broadening rejection correction
  cxs // Constant cross section
};

// Electron treatments
// TODO: Convert to enum
constexpr int ELECTRON_LED {1}; // Local Energy Deposition
constexpr int ELECTRON_TTB {2}; // Thick Target Bremsstrahlung

// ============================================================================
// MULTIGROUP RELATED

// MGXS Table Types
// TODO: Convert to enum
constexpr int MGXS_ISOTROPIC {1}; // Isotroically weighted data
constexpr int MGXS_ANGLE {2};     // Data by angular bins

// Flag to denote this was a macroscopic data object
constexpr double MACROSCOPIC_AWR {-2.};

// Number of mu bins to use when converting Legendres to tabular type
constexpr int DEFAULT_NMU {33};

// Mgxs::get_xs enumerated types
// TODO: Convert to enum
constexpr int MG_GET_XS_TOTAL              {0};
constexpr int MG_GET_XS_ABSORPTION         {1};
constexpr int MG_GET_XS_INVERSE_VELOCITY   {2};
constexpr int MG_GET_XS_DECAY_RATE         {3};
constexpr int MG_GET_XS_SCATTER            {4};
constexpr int MG_GET_XS_SCATTER_MULT       {5};
constexpr int MG_GET_XS_SCATTER_FMU_MULT   {6};
constexpr int MG_GET_XS_SCATTER_FMU        {7};
constexpr int MG_GET_XS_FISSION            {8};
constexpr int MG_GET_XS_KAPPA_FISSION      {9};
constexpr int MG_GET_XS_PROMPT_NU_FISSION  {10};
constexpr int MG_GET_XS_DELAYED_NU_FISSION {11};
constexpr int MG_GET_XS_NU_FISSION         {12};
constexpr int MG_GET_XS_CHI_PROMPT         {13};
constexpr int MG_GET_XS_CHI_DELAYED        {14};

// ============================================================================
// TALLY-RELATED CONSTANTS

// Tally result entries
constexpr int RESULT_VALUE  {0};
constexpr int RESULT_SUM    {1};
constexpr int RESULT_SUM_SQ {2};

// Tally type
// TODO: Convert to enum
constexpr int TALLY_VOLUME          {1};
constexpr int TALLY_MESH_SURFACE    {2};
constexpr int TALLY_SURFACE         {3};

// Tally estimator types
// TODO: Convert to enum
constexpr int ESTIMATOR_ANALOG      {1};
constexpr int ESTIMATOR_TRACKLENGTH {2};
constexpr int ESTIMATOR_COLLISION   {3};

// Event types for tallies
// TODO: Convert to enum
constexpr int EVENT_SURFACE {-2};
constexpr int EVENT_LATTICE {-1};
constexpr int EVENT_KILL    {0};
constexpr int EVENT_SCATTER {1};
constexpr int EVENT_ABSORB  {2};

// Tally score type -- if you change these, make sure you also update the
// _SCORES dictionary in openmc/capi/tally.py
// TODO: Convert to enum
constexpr int SCORE_FLUX               {-1}; // flux
constexpr int SCORE_TOTAL              {-2}; // total reaction rate
constexpr int SCORE_SCATTER            {-3}; // scattering rate
constexpr int SCORE_NU_SCATTER         {-4}; // scattering production rate
constexpr int SCORE_ABSORPTION         {-5}; // absorption rate
constexpr int SCORE_FISSION            {-6}; // fission rate
constexpr int SCORE_NU_FISSION         {-7}; // neutron production rate
constexpr int SCORE_KAPPA_FISSION      {-8}; // fission energy production rate
constexpr int SCORE_CURRENT            {-9}; // current
constexpr int SCORE_EVENTS             {-10}; // number of events
constexpr int SCORE_DELAYED_NU_FISSION {-11}; // delayed neutron production rate
constexpr int SCORE_PROMPT_NU_FISSION  {-12}; // prompt neutron production rate
constexpr int SCORE_INVERSE_VELOCITY   {-13}; // flux-weighted inverse velocity
constexpr int SCORE_FISS_Q_PROMPT      {-14}; // prompt fission Q-value
constexpr int SCORE_FISS_Q_RECOV       {-15}; // recoverable fission Q-value
constexpr int SCORE_DECAY_RATE         {-16}; // delayed neutron precursor decay rate
constexpr int SCORE_HEATING            {-17}; // nuclear heating (neutron or photon)

// Tally map bin finding
constexpr int NO_BIN_FOUND {-1};

// Tally filter and map types
// TODO: Refactor to remove or convert to enum
constexpr int FILTER_UNIVERSE       {1};
constexpr int FILTER_MATERIAL       {2};
constexpr int FILTER_CELL           {3};

// Mesh types
constexpr int MESH_REGULAR {1};

// Tally surface current directions
constexpr int OUT_LEFT   {1};  // x min
constexpr int IN_LEFT    {2};  // x min
constexpr int OUT_RIGHT  {3};  // x max
constexpr int IN_RIGHT   {4};  // x max
constexpr int OUT_BACK   {5};  // y min
constexpr int IN_BACK    {6};  // y min
constexpr int OUT_FRONT  {7};  // y max
constexpr int IN_FRONT   {8};  // y max
constexpr int OUT_BOTTOM {9};  // z min
constexpr int IN_BOTTOM  {10}; // z min
constexpr int OUT_TOP    {11}; // z max
constexpr int IN_TOP     {12}; // z max

// Global tally parameters
constexpr int N_GLOBAL_TALLIES {4};
constexpr int K_COLLISION   {0};
constexpr int K_ABSORPTION  {1};
constexpr int K_TRACKLENGTH {2};
constexpr int LEAKAGE       {3};

// Miscellaneous
constexpr int C_NONE {-1};
constexpr int F90_NONE {0}; //TODO: replace usage of this with C_NONE

// Interpolation rules
enum class Interpolation {
  histogram = 1, lin_lin = 2, lin_log = 3, log_lin = 4, log_log = 5
};

// Run modes
constexpr int RUN_MODE_FIXEDSOURCE {1};
constexpr int RUN_MODE_EIGENVALUE  {2};
constexpr int RUN_MODE_PLOTTING    {3};
constexpr int RUN_MODE_PARTICLE    {4};
constexpr int RUN_MODE_VOLUME      {5};

// ============================================================================
// CMFD CONSTANTS

// For non-accelerated regions on coarse mesh overlay
constexpr int CMFD_NOACCEL {-1};

} // namespace openmc

#endif // OPENMC_CONSTANTS_H
