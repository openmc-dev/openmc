#ifndef OPENMC_TALLIES_FILTER_UNIVERSE_H
#define OPENMC_TALLIES_FILTER_UNIVERSE_H

#include <cstdint>
#include <unordered_map>
#include <vector>

#include <gsl/gsl>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Specifies which geometric universes tally events reside in.
//==============================================================================

class UniverseFilter : public Filter
{
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~UniverseFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "universe";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  void set_universes(gsl::span<int32_t> universes);

private:
  //----------------------------------------------------------------------------
  // Data members

  //! The indices of the universes binned by this filter.
  std::vector<int32_t> universes_;

  //! A map from universe indices to filter bin indices.
  std::unordered_map<int32_t, int> map_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_UNIVERSE_H
