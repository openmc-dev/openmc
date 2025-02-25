#ifndef OPENMC_NCRYSTAL_INTERFACE_H
#define OPENMC_NCRYSTAL_INTERFACE_H

#include "openmc/ncrystal_load.h"
#include "openmc/particle.h"

#include <cstdint> // for uint64_t
#include <limits>  // for numeric_limits
#include <string>

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

//! Energy in [eV] to switch between NCrystal and ENDF
constexpr double NCRYSTAL_MAX_ENERGY {5.0};

//==============================================================================
// Wrapper class an NCrystal material
//==============================================================================

class NCrystalMat {
public:
  //----------------------------------------------------------------------------
  // Constructors
  NCrystalMat() = default;//empty object
  explicit NCrystalMat(const std::string& cfg);

  //----------------------------------------------------------------------------
  // Methods

  //! Return configuration string:
  const std::string& cfg() const { return cfg_; }

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
  operator bool() const { return !cfg_.empty(); }

  NCrystalMat clone() const
  {
    NCrystalMat c;
    c.cfg_ = cfg_;
    c.proc_ = proc_.clone();
    return c;
  }

private:
  //----------------------------------------------------------------------------
  // Data members (only present when compiling with NCrystal support)
  std::string cfg_; //!< NCrystal configuration string
  NCrystalScatProc proc_; //!< NCrystal scatter process
};

//==============================================================================
// Functions
//==============================================================================

void ncrystal_update_micro(double xs, NuclideMicroXS& micro);

} // namespace openmc

#endif // OPENMC_NCRYSTAL_INTERFACE_H
