#ifndef OPENMC_TALLIES_FILTER_TIME_H
#define OPENMC_TALLIES_FILTER_TIME_H

#include "openmc/span.h"
#include "openmc/tallies/filter.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Bins the incident particle time.
//==============================================================================

class TimeFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~TimeFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "time"; }
  FilterType type() const override { return FilterType::TIME; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const vector<double>& bins() const { return bins_; }
  void set_bins(span<const double> bins);

protected:
  //----------------------------------------------------------------------------
  // Data members

  vector<double> bins_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_ENERGY_H
