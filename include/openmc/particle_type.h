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
#include "openmc/error.h"

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

  // Constructor from PDG number
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

  // Get particle mass in [u]
  double mass() const
  {
    int32_t p = std::abs(pdg_number_);
    if (ATOMIC_MASS.count(p)) {
      return ATOMIC_MASS[p];
    } else {
      fatal_error("Unknown mass for particle " + str());
    }
  }

  // Convert to string representation
  std::string str() const;

  // Check if this represents a nucleus (vs elementary particle)
  constexpr bool is_nucleus() const
  {
    // PDG nuclear codes are >= 1000000000 (100ZZZAAAI format)
    return pdg_number_ >= 1000000000;
  }

  // Get transport index (0-3 for transportable particles, C_NONE otherwise)
  constexpr int transport_index() const;

  // Check if this is a neutron
  constexpr bool is_neutron() const { return pdg_number_ == PDG_NEUTRON; }

  // Check if this is a photon
  constexpr bool is_photon() const { return pdg_number_ == PDG_PHOTON; }

  constexpr bool is_transportable() const
  {
    return this->transport_index() != C_NONE;
  }

  //----------------------------------------------------------------------------
  // Static factory methods

  static constexpr ParticleType neutron() { return ParticleType {PDG_NEUTRON}; }
  static constexpr ParticleType photon() { return ParticleType {PDG_PHOTON}; }
  static constexpr ParticleType electron()
  {
    return ParticleType {PDG_ELECTRON};
  }
  static constexpr ParticleType positron()
  {
    return ParticleType {PDG_POSITRON};
  }
  static constexpr ParticleType proton() { return ParticleType {PDG_PROTON}; }
  static constexpr ParticleType deuteron()
  {
    return ParticleType {PDG_DEUTERON};
  }
  static constexpr ParticleType triton() { return ParticleType {PDG_TRITON}; }
  static constexpr ParticleType alpha() { return ParticleType {PDG_ALPHA}; }

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
// ParticleType member function implementations (inline)
//------------------------------------------------------------------------------

constexpr int ParticleType::transport_index() const
{
  switch (pdg_number_) {
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

//------------------------------------------------------------------------------
// Legacy conversion helpers
//------------------------------------------------------------------------------

// Legacy enum code (0..3) to ParticleType conversion
ParticleType legacy_particle_index_to_type(int code);

} // namespace openmc

#endif // OPENMC_PARTICLE_TYPE_H
