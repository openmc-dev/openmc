#ifndef OPENMC_TALLIES_FILTER_LIFETIME_MOMENT_H
#define OPENMC_TALLIES_FILTER_LIFETIME_MOMENT_H

#include <string>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Gives power moments of the time elapsed since the particle was born
//!
//! The weight of bin n (n = 0, ..., order) is tau^n, where tau is the time in
//! seconds since the particle was started from its source site (see
//! Particle::lifetime). Secondary particles of (n,xn) reactions and split
//! particles restart tau at the collision or split that created them, a
//! limitation shared with MeshBornFilter and the IFP lifetime. Because the
//! weight varies along a track, tallies using this filter are scored with a
//! collision or analog estimator.
//==============================================================================

class LifetimeMomentFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  LifetimeMomentFilter() { this->set_order(0); }

  ~LifetimeMomentFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "lifetimemoment"; }
  FilterType type() const override { return FilterType::LIFETIME_MOMENT; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  int order() const { return order_; }

  void set_order(int order);

private:
  //----------------------------------------------------------------------------
  // Data members

  int order_ {0}; //!< Highest moment order N; bins are n = 0, ..., N
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_LIFETIME_MOMENT_H
