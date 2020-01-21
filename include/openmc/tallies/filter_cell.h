#ifndef OPENMC_TALLIES_FILTER_CELL_H
#define OPENMC_TALLIES_FILTER_CELL_H

#include <cstdint>
#include <unordered_map>
#include <vector>

#include <gsl/gsl>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Specifies which geometric cells tally events reside in.
//==============================================================================

class CellFilter : public Filter
{
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~CellFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "cell";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const std::vector<int32_t>& cells() const { return cells_; }

  void set_cells(gsl::span<int32_t> cells);

protected:
  //----------------------------------------------------------------------------
  // Data members

  //! The indices of the cells binned by this filter.
  std::vector<int32_t> cells_;

  //! A map from cell indices to filter bin indices.
  std::unordered_map<int32_t, int> map_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_CELL_H
