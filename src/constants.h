//! \file constants.h
//! A collection of constants

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <array>
#include <vector>

namespace openmc {

typedef std::array<double, 3> dir_arr;
typedef std::vector<double> double_1dvec;
typedef std::vector<std::vector<double> > double_2dvec;
typedef std::vector<std::vector<std::vector<double> > > double_3dvec;
typedef std::vector<std::vector<std::vector<std::vector<double> > > > double_4dvec;
typedef std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > double_5dvec;
typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > > double_6dvec;
typedef std::vector<int> int_1dvec;
typedef std::vector<std::vector<int> > int_2dvec;
typedef std::vector<std::vector<std::vector<int> > > int_3dvec;

int constexpr MAX_SAMPLE {10000};

constexpr std::array<int, 3> VERSION {0, 10, 0};
constexpr std::array<int, 2> VERSION_PARTICLE_RESTART {2, 0};

// Maximum number of words in a single line, length of line, and length of
// single word
constexpr int MAX_WORDS {500};
constexpr int MAX_LINE_LEN {250};
constexpr int MAX_WORD_LEN {150};
constexpr int MAX_FILE_LEN {255};

// Physical Constants
constexpr double K_BOLTZMANN {8.6173303e-5};  // Boltzmann constant in eV/K

// Angular distribution type
constexpr int ANGLE_ISOTROPIC {1};
constexpr int ANGLE_32_EQUI {2};
constexpr int ANGLE_TABULAR {3};
constexpr int ANGLE_LEGENDRE {4};
constexpr int ANGLE_HISTOGRAM {5};

// MGXS Table Types
constexpr int MGXS_ISOTROPIC {1}; // Isotroically weighted data
constexpr int MGXS_ANGLE {2};     // Data by angular bins

// Flag to denote this was a macroscopic data object
constexpr double MACROSCOPIC_AWR {-2.};

// Number of mu bins to use when converting Legendres to tabular type
constexpr int DEFAULT_NMU {33};

// Temperature treatment method
constexpr int TEMPERATURE_NEAREST {1};
constexpr int TEMPERATURE_INTERPOLATION {2};

// TODO: cmath::M_PI has 3 more digits precision than the Fortran constant we
// use so for now we will reuse the Fortran constant until we are OK with
// modifying test results
constexpr double PI {3.1415926535898};

const double SQRT_PI {std::sqrt(PI)};


} // namespace openmc

#endif // CONSTANTS_H