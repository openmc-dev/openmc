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

  //! Record, for each of this filter's bins, which entry of
  //! model::active_point_detectors it sits on. Called once per batch from
  //! setup_active_tallies(), after that list has been assembled.
  void build_detector_bins();

private:
  //----------------------------------------------------------------------------
  // Data members

  vector<std::pair<Position, double>> detectors_;

  //! Parallel to detectors_: the index into model::active_point_detectors
  //! that each bin sits on, or C_NONE if this filter's tally is not active.
  //! Mapping this way rather than the reverse keeps the behaviour of a filter
  //! that puts two bins on one position -- different exclusion radii at the
  //! same point -- where a detector-to-bin map could only name one of them.
  //! Rebuilt every batch by build_detector_bins().
  vector<int> bin_detector_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_POINT_H
