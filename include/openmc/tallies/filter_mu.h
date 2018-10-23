#ifndef OPENMC_TALLIES_FILTER_MU_H
#define OPENMC_TALLIES_FILTER_MU_H

#include <vector>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Bins the incoming-outgoing direction cosine.  This is only used for scatter
//! reactions.
//==============================================================================

class MuFilter : public Filter
{
public:
  ~MuFilter() = default;

  std::string type() const override {return "mu";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  std::vector<double> bins_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_MU_H
