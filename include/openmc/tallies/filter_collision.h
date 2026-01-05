#ifndef OPENMC_TALLIES_FILTER_COLLISIONS_H
#define OPENMC_TALLIES_FILTER_COLLISIONS_H

#include <unordered_map>

#include "openmc/span.h"
#include "openmc/tallies/filter.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Bins the incident neutron energy.
//==============================================================================

class CollisionFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~CollisionFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "collision"; }
  FilterType type() const override { return FilterType::COLLISION; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const vector<int>& bins() const { return bins_; }
  void set_bins(span<const int> bins);

protected:
  //----------------------------------------------------------------------------
  // Data members

  vector<int> bins_;

  std::unordered_map<int, int> map_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_COLLISIONS_H
