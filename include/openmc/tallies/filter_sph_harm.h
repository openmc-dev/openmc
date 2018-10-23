#ifndef OPENMC_TALLIES_FILTER_SPH_HARM_H
#define OPENMC_TALLIES_FILTER_SPH_HARM_H

#include <string>

#include "openmc/tallies/filter.h"

namespace openmc {

//TODO: those integer values are not needed when Fortran interop is removed
enum class SphericalHarmonicsCosine {
  scatter = 1, particle = 2
};

//==============================================================================
//! Gives spherical harmonics expansion moments of a tally score
//==============================================================================

class SphericalHarmonicsFilter : public Filter
{
public:
  ~SphericalHarmonicsFilter() = default;

  std::string type() const override {return "sphericalharmonics";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  int order_;

  //! The type of angle that this filter measures when binning events.
  SphericalHarmonicsCosine cosine_ {SphericalHarmonicsCosine::particle};
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_SPH_HARM_H
