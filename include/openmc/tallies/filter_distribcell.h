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
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~DistribcellFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "distribcell";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  int32_t cell() const { return cell_; }

  void set_cell(int32_t cell);

private:
  //----------------------------------------------------------------------------
  // Data members

  int32_t cell_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_DISTRIBCELL_H
