//==============================================================================
// atomic masses definitions
//==============================================================================

#ifndef OPENMC_ATOMIC_MASS_H
#define OPENMC_ATOMIC_MASS_H

#include <cstdint>
#include <unordered_map>

namespace openmc {

// The bare-particle rest masses below are from the Committee on Data for
// Science and Technology (CODATA) 2018 recommendation
// (https://physics.nist.gov/cuu/Constants/). The ATOMIC_MASS table declared
// further down holds neutral atomic masses from AME2020 instead.

// Physical constants
constexpr double MASS_ELECTRON {5.48579909065e-4}; // mass of an electron in amu
constexpr double MASS_NEUTRON {1.00866491595};     // mass of a neutron in amu
constexpr double MASS_PROTON {1.007276466621};     // mass of a proton in amu
constexpr double MASS_DEUTRON {2.013553212745};    // mass of a deuteron in amu
constexpr double MASS_TRITON {3.01550071621};      // mass of a triton in amu
constexpr double MASS_HELION {3.014932247175};     // mass of a helion in amu
constexpr double MASS_ALPHA {4.001506179127};      // mass of an alpha in amu

//! Neutral ground-state atomic masses in [u], indexed by nuclear PDG code
extern const std::unordered_map<int32_t, double> ATOMIC_MASS;

//! Return the neutral ground-state atomic mass for a PDG code in [u]
//!
//! Nuclear isomer codes are normalized to the ground state, and the proton
//! code is treated as an alias for H-1. Returns zero if the code does not
//! identify a nuclide or its mass is not tabulated.
double atomic_mass_from_pdg(int32_t pdg);

//! Return the neutral ground-state atomic mass of a nuclide in [u]
//!
//! Returns zero if the nuclide is invalid or its mass is not tabulated.
double atomic_mass(int Z, int A);

//! Return the bare-particle rest mass for a PDG code in [u]
//!
//! Elementary particles and light nuclei use their CODATA masses. Other
//! nuclei are derived from the corresponding neutral atomic mass. Returns zero
//! for a photon or if the particle mass is not available.
double nuclear_mass_from_pdg(int32_t pdg);

//! Return the nuclear mass of a nuclide in [u]
//!
//! For light particles, the CODATA bare-particle mass is returned. Otherwise,
//! the mass is approximated by subtracting the masses of the atomic electrons
//! from the tabulated atomic mass. Electron binding energy is neglected. A
//! nuclide with Z = 0 and A = 1 is the free neutron, for which the CODATA
//! neutron mass is returned; atomic_mass() reports zero for it instead, as
//! there is no corresponding neutral atom.
//! Returns zero if the nuclide is invalid or its atomic mass is not tabulated.
double nuclear_mass(int Z, int A);

} // namespace openmc

#endif //  OPENMC_ATOMIC_MASS_H
