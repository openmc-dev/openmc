#ifndef OPENMC_TALLIES_FILTER_SPH_HARM_H
#define OPENMC_TALLIES_FILTER_SPH_HARM_H

#include <string>

#include <gsl/gsl>

#include "openmc/tallies/filter.h"

namespace openmc {

enum class SphericalHarmonicsCosine {
  scatter, particle
};

//==============================================================================
//! Gives spherical harmonics expansion moments of a tally score
//==============================================================================

class SphericalHarmonicsFilter : public Filter
{
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~SphericalHarmonicsFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "sphericalharmonics";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  int order() const { return order_; }

  void set_order(int order);

  SphericalHarmonicsCosine cosine() const { return cosine_; }

  void set_cosine(gsl::cstring_span cosine);

private:
  //----------------------------------------------------------------------------
  // Data members

  int order_;

  //! The type of angle that this filter measures when binning events.
  SphericalHarmonicsCosine cosine_ {SphericalHarmonicsCosine::particle};
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_SPH_HARM_H
