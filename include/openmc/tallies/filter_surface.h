#ifndef OPENMC_TALLIES_FILTER_SURFACE_H
#define OPENMC_TALLIES_FILTER_SURFACE_H

#include <cstdint>
#include <unordered_map>
#include <vector>

#include <gsl/gsl>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Specifies which surface particles are crossing
//==============================================================================

class SurfaceFilter : public Filter
{
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~SurfaceFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "surface";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  void set_surfaces(gsl::span<int32_t> surfaces);

private:
  //----------------------------------------------------------------------------
  // Data members

  //! The indices of the surfaces binned by this filter.
  std::vector<int32_t> surfaces_;

  //! A map from surface indices to filter bin indices.
  std::unordered_map<int32_t, int> map_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_SURFACE_H
