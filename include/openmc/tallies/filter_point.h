#ifndef OPENMC_TALLIES_FILTER_POINT_H
#define OPENMC_TALLIES_FILTER_POINT_H

#include "openmc/position.h"
#include "openmc/span.h"
#include "openmc/tallies/filter.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Bins tally by point detectors
//==============================================================================

class PointFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~PointFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "point"; }
  FilterType type() const override { return FilterType::POINT; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const vector<std::pair<Position, double>>& detectors() const
  {
    return detectors_;
  }

  void set_detectors(span<std::pair<Position, double>> detectors);

private:
  //----------------------------------------------------------------------------
  // Data members

  vector<std::pair<Position, double>> detectors_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_PARTICLE_H
