//==============================================================================
// atomic masses definitions
//==============================================================================

#ifndef OPENMC_ATOMIC_MASS_H
#define OPENMC_ATOMIC_MASS_H

#include <cstdint>
#include <unordered_map>

namespace openmc {

// Values here are from the Committee on Data for Science and Technology
// (CODATA) 2018 recommendation (https://physics.nist.gov/cuu/Constants/).

// Physical constants
constexpr double MASS_ELECTRON {5.48579909065e-4}; // mass of an electron in amu
constexpr double MASS_NEUTRON {1.00866491595};     // mass of a neutron in amu
constexpr double MASS_PROTON {1.007276466621};     // mass of a proton in amu
constexpr double MASS_DEUTRON {2.013553212745};    // mass of a deutron in amu
constexpr double MASS_HELION {3.014932247175};     // mass of a helion in amu
constexpr double MASS_ALPHA {4.001506179127};      // mass of an alpha in amu

extern std::unordered_map<int32_t, double> ATOMIC_MASS;

} // namespace openmc

#endif //  OPENMC_ATOMIC_MASS_H
