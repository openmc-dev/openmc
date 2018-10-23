#ifndef OPENMC_TALLIES_FILTER_CELL_H
#define OPENMC_TALLIES_FILTER_CELL_H

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Specifies which geometric cells tally events reside in.
//==============================================================================

class CellFilter : public Filter
{
public:
  ~CellFilter() = default;

  std::string type() const override {return "cell";}

  void from_xml(pugi::xml_node node) override;

  void initialize() override;

  void get_all_bins(Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  std::vector<int32_t> cells_;
  std::unordered_map<int32_t, int> map_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_CELL_H
