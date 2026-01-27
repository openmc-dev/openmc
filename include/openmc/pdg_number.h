//==============================================================================
// Particle PDG definitions and helpers
//==============================================================================

#ifndef OPENMC_PDG_NUMBER_H
#define OPENMC_PDG_NUMBER_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <type_traits>

#include "openmc/constants.h"

namespace openmc {

//------------------------------------------------------------------------------
// Canonical PDG particle identifier (standard-layout, trivially copyable)
//------------------------------------------------------------------------------

struct PDGNumber {
  int32_t value;
};

static_assert(
  std::is_standard_layout_v<PDGNumber>, "ParticlePdg must be standard-layout");
static_assert(std::is_trivially_copyable_v<PDGNumber>,
  "ParticlePdg must be trivially copyable");
static_assert(
  offsetof(PDGNumber, value) == 0, "ParticlePdg value must be at offset 0");
static_assert(sizeof(PDGNumber) == sizeof(int32_t),
  "ParticlePdg must be same size as int32_t");

constexpr bool operator==(PDGNumber lhs, PDGNumber rhs)
{
  return lhs.value == rhs.value;
}
constexpr bool operator!=(PDGNumber lhs, PDGNumber rhs)
{
  return lhs.value != rhs.value;
}
constexpr bool operator<(PDGNumber lhs, PDGNumber rhs)
{
  return lhs.value < rhs.value;
}

//------------------------------------------------------------------------------
// PDG constants (canonical particle identity)
//------------------------------------------------------------------------------

inline constexpr PDGNumber PDG_NEUTRON {2112};
inline constexpr PDGNumber PDG_PHOTON {22};
inline constexpr PDGNumber PDG_ELECTRON {11};
inline constexpr PDGNumber PDG_POSITRON {-11};
inline constexpr PDGNumber PDG_PROTON {2212};
inline constexpr PDGNumber PDG_DEUTERON {1000010020};
inline constexpr PDGNumber PDG_TRITON {1000010030};
inline constexpr PDGNumber PDG_ALPHA {1000020040};

//------------------------------------------------------------------------------
// Transport indexing helpers (for dense arrays keyed by transport kind)
//------------------------------------------------------------------------------

constexpr int transport_index_from_pdg(PDGNumber pdg)
{
  switch (pdg.value) {
  case PDG_NEUTRON.value:
    return 0;
  case PDG_PHOTON.value:
    return 1;
  case PDG_ELECTRON.value:
    return 2;
  case PDG_POSITRON.value:
    return 3;
  default:
    return C_NONE;
  }
}

constexpr bool is_transport_pdg(PDGNumber pdg)
{
  return transport_index_from_pdg(pdg) != C_NONE;
}

//------------------------------------------------------------------------------
// Conversion helpers
//------------------------------------------------------------------------------

// Convert a particle identifier string to PDG.
PDGNumber str_to_pdg_number(std::string_view str);

// Convert PDG to a canonical string (known names or GNDS nuclide name).
std::string pdg_number_to_str(PDGNumber pdg);

// Legacy enum code (0..3) to PDG conversion
PDGNumber legacy_particle_code_to_pdg(int code);

// GNDS nuclide name to PDG nuclear code (100ZZZAAAI)
PDGNumber pdg_number_from_nuclide_name(std::string_view gnds);

// PDG nuclear code to GNDS nuclide name (e.g., "Fe57", "Am242_m1")
std::string nuclide_name_from_pdg_number(PDGNumber pdg);

// Build a PDG nuclear code from Z, A, and isomer level
PDGNumber pdg_number_from_zam(int Z, int A, int m);

// Check for PDG nuclear code (100ZZZAAAI)
bool is_nuclear_pdg(PDGNumber pdg);

} // namespace openmc

#endif // OPENMC_PDG_NUMBER_H
