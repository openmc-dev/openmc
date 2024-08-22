#include "openmc/ncrystal_interface.h"

#include "openmc/error.h"
#include "openmc/material.h"
#include "openmc/random_lcg.h"

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

#ifdef NCRYSTAL
const bool OPENMC_API NCRYSTAL_ENABLED = true;
#else
const bool OPENMC_API NCRYSTAL_ENABLED = false;
#endif

//==============================================================================
// NCrystal wrapper class for the OpenMC random number generator
//==============================================================================

#ifdef NCRYSTAL
class NCrystalRNGWrapper : public NCrystal::RNGStream {
public:
  constexpr NCrystalRNGWrapper(uint64_t* seed) noexcept : openmc_seed_(seed) {}

protected:
  double actualGenerate() override
  {
    return std::max<double>(
      std::numeric_limits<double>::min(), prn(openmc_seed_));
  }

private:
  uint64_t* openmc_seed_;
};
#endif

//==============================================================================
// NCrystal implementation
//==============================================================================

NCrystalMat::NCrystalMat(const std::string& cfg)
{
#ifdef NCRYSTAL
  cfg_ = cfg;
  ptr_ = NCrystal::FactImpl::createScatter(cfg);
#else
  fatal_error("Your build of OpenMC does not support NCrystal materials.");
#endif
}

#ifdef NCRYSTAL
std::string NCrystalMat::cfg() const
{
  return cfg_;
}

double NCrystalMat::xs(const Particle& p) const
{
  // Calculate scattering XS per atom with NCrystal, only once per material
  NCrystal::CachePtr dummy_cache;
  auto nc_energy = NCrystal::NeutronEnergy {p.E()};
  return ptr_->crossSection(dummy_cache, nc_energy, {p.u().x, p.u().y, p.u().z})
    .get();
}

void NCrystalMat::scatter(Particle& p) const
{
  NCrystalRNGWrapper rng(p.current_seed()); // Initialize RNG
  // create a cache pointer for multi thread physics
  NCrystal::CachePtr dummy_cache;
  auto nc_energy = NCrystal::NeutronEnergy {p.E()};
  auto outcome = ptr_->sampleScatter(
    dummy_cache, rng, nc_energy, {p.u().x, p.u().y, p.u().z});

  // Modify attributes of particle
  p.E() = outcome.ekin.get();
  Direction u_old {p.u()};
  p.u() =
    Direction(outcome.direction[0], outcome.direction[1], outcome.direction[2]);
  p.mu() = u_old.dot(p.u());
  p.event_mt() = ELASTIC;
}

NCrystalMat::operator bool() const
{
  return ptr_.get();
}
#endif

//==============================================================================
// Functions
//==============================================================================

void ncrystal_update_micro(double xs, NuclideMicroXS& micro)
{
  if (micro.thermal > 0 || micro.thermal_elastic > 0) {
    fatal_error("S(a,b) treatment and NCrystal are not compatible.");
  }
  // remove free atom cross section
  // and replace it by scattering cross section per atom from NCrystal
  micro.total = micro.total - micro.elastic + xs;
  micro.elastic = xs;
}

} // namespace openmc
