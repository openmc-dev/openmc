#ifndef OPENMC_TALLIES_FILTER_DISTRIBCELL_H
#define OPENMC_TALLIES_FILTER_DISTRIBCELL_H

#include <string>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Specifies which distributed geometric cells tally events reside in.
//==============================================================================

class DistribcellFilter : public Filter
{
public:
  ~DistribcellFilter() = default;

  std::string type() const override {return "distribcell";}

  void from_xml(pugi::xml_node node) override;

  void initialize() override;

  void get_all_bins(Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  int32_t cell_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_DISTRIBCELL_H
