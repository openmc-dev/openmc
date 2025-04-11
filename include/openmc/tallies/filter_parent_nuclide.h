#ifndef OPENMC_TALLIES_FILTER_PARENT_NUCLIDE_H
#define OPENMC_TALLIES_FILTER_PARENT_NUCLIDE_H

#include <string>
#include <unordered_map>

#include "openmc/span.h"
#include "openmc/tallies/filter.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Bins events by parent nuclide (for decay photons)
//==============================================================================

class ParentNuclideFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~ParentNuclideFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "parentnuclide"; }
  FilterType type() const override { return FilterType::PARENT_NUCLIDE; }

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
  vector<std::string> nuclides_;

  std::unordered_map<int, int> map_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_PARENT_NUCLIDE_H
