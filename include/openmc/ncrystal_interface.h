#ifndef OPENMC_NCRYSTAL_INTERFACE_H
#define OPENMC_NCRYSTAL_INTERFACE_H

#ifdef NCRYSTAL
#include "NCrystal/NCrystal.hh"
#endif

#include "openmc/particle.h"

#include <cstdint> // for uint64_t
#include <limits>  // for numeric_limits
#include <string>

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

extern "C" const bool NCRYSTAL_ENABLED;

//! Energy in [eV] to switch between NCrystal and ENDF
constexpr double NCRYSTAL_MAX_ENERGY {5.0};

//==============================================================================
// Wrapper class an NCrystal material
//==============================================================================

class NCrystalMat {
public:
  //----------------------------------------------------------------------------
  // Constructors
  NCrystalMat() = default;
  explicit NCrystalMat(const std::string& cfg);

  //----------------------------------------------------------------------------
  // Methods

#ifdef NCRYSTAL
  //! Return configuration string
  std::string cfg() const;

  //! Get cross section from NCrystal material
  //
  //! \param[in] p  Particle object
  //! \return  Cross section in [b]
  double xs(const Particle& p) const;

  // Process scattering event
  //
  //! \param[in] p  Particle object
  void scatter(Particle& p) const;

  //! Whether the object holds a valid NCrystal material
  operator bool() const;
#else

  //----------------------------------------------------------------------------
  // Trivial methods when compiling without NCRYSTAL
  std::string cfg() const
  {
    return "";
  }
  double xs(const Particle& p) const
  {
    return -1.0;
  }
  void scatter(Particle& p) const {}
  operator bool() const
  {
    return false;
  }
#endif

private:
  //----------------------------------------------------------------------------
  // Data members (only present when compiling with NCrystal support)
#ifdef NCRYSTAL
  std::string cfg_; //!< NCrystal configuration string
  std::shared_ptr<const NCrystal::ProcImpl::Process>
    ptr_; //!< Pointer to NCrystal material object
#endif
};

//==============================================================================
// Functions
//==============================================================================

void ncrystal_update_micro(double xs, NuclideMicroXS& micro);

} // namespace openmc

#endif // OPENMC_NCRYSTAL_INTERFACE_H
