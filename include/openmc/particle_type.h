//==============================================================================
// ParticleType class definition
//==============================================================================

#ifndef OPENMC_PARTICLE_TYPE_H
#define OPENMC_PARTICLE_TYPE_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <type_traits>

#include "openmc/constants.h"

namespace openmc {

//------------------------------------------------------------------------------
// PDG constants (canonical particle identity as simple integers)
//------------------------------------------------------------------------------

inline constexpr int32_t PDG_NEUTRON = 2112;
inline constexpr int32_t PDG_PHOTON = 22;
inline constexpr int32_t PDG_ELECTRON = 11;
inline constexpr int32_t PDG_POSITRON = -11;
inline constexpr int32_t PDG_PROTON = 2212;
inline constexpr int32_t PDG_DEUTERON = 1000010020;
inline constexpr int32_t PDG_TRITON = 1000010030;
inline constexpr int32_t PDG_ALPHA = 1000020040;

//------------------------------------------------------------------------------
// ParticleType class (standard-layout, trivially copyable)
//------------------------------------------------------------------------------

class ParticleType {
public:
  //----------------------------------------------------------------------------
  // Constructors

  // Default constructor: defaults to neutron
  constexpr ParticleType() : pdg_number_(PDG_NEUTRON) {}

  // Constructor from PDG code
  constexpr explicit ParticleType(int32_t pdg_number) : pdg_number_(pdg_number)
  {}

  // Constructor from particle name string (e.g., "neutron", "photon", "Fe56")
  explicit ParticleType(std::string_view str);

  // Constructor from Z, A, and metastable state for nuclear particles
  constexpr ParticleType(int Z, int A, int m = 0)
    : pdg_number_(1000000000 + Z * 10000 + A * 10 + m)
  {}

  //----------------------------------------------------------------------------
  // Accessors

  // Accessor for the underlying PDG number
  constexpr int32_t pdg_number() const { return pdg_number_; }

  //----------------------------------------------------------------------------
  // Methods

  // Convert to string representation
  std::string str() const;

  // Check if this represents a nucleus (vs elementary particle)
  constexpr bool is_nucleus() const
  {
    // PDG nuclear codes are >= 1000000000 (100ZZZAAAI format)
    return pdg_number_ >= 1000000000;
  }

private:
  int32_t pdg_number_;
};

//------------------------------------------------------------------------------
// Static assertions to ensure standard-layout and trivially copyable
//------------------------------------------------------------------------------

static_assert(std::is_standard_layout_v<ParticleType>,
  "ParticleType must be standard-layout");
static_assert(std::is_trivially_copyable_v<ParticleType>,
  "ParticleType must be trivially copyable");
static_assert(sizeof(ParticleType) == sizeof(int32_t),
  "ParticleType must be same size as int32_t");

//------------------------------------------------------------------------------
// Comparison operators (free functions for symmetry)
//------------------------------------------------------------------------------

constexpr bool operator==(ParticleType lhs, ParticleType rhs)
{
  return lhs.pdg_number() == rhs.pdg_number();
}

constexpr bool operator!=(ParticleType lhs, ParticleType rhs)
{
  return lhs.pdg_number() != rhs.pdg_number();
}

constexpr bool operator<(ParticleType lhs, ParticleType rhs)
{
  return lhs.pdg_number() < rhs.pdg_number();
}

//------------------------------------------------------------------------------
// Transport indexing helpers (for dense arrays keyed by transport kind)
//------------------------------------------------------------------------------

constexpr int transport_index(ParticleType type)
{
  switch (type.pdg_number()) {
  case PDG_NEUTRON:
    return 0;
  case PDG_PHOTON:
    return 1;
  case PDG_ELECTRON:
    return 2;
  case PDG_POSITRON:
    return 3;
  default:
    return C_NONE;
  }
}

constexpr bool is_transportable(ParticleType type)
{
  return transport_index(type) != C_NONE;
}

//------------------------------------------------------------------------------
// Legacy conversion helpers
//------------------------------------------------------------------------------

// Legacy enum code (0..3) to ParticleType conversion
ParticleType legacy_particle_index_to_type(int code);

} // namespace openmc

#endif // OPENMC_PARTICLE_TYPE_H
