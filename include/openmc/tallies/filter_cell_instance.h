#ifndef OPENMC_TALLIES_FILTER_CELL_INSTANCE_H
#define OPENMC_TALLIES_FILTER_CELL_INSTANCE_H

#include <cstdint>
#include <unordered_map>

#include <gsl/gsl-lite.hpp>

#include "openmc/cell.h"
#include "openmc/tallies/filter.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Specifies cell instances that tally events reside in.
//==============================================================================

class CellInstanceFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  CellInstanceFilter() = default;
  CellInstanceFilter(gsl::span<CellInstance> instances);
  ~CellInstanceFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override { return "cellinstance"; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const vector<CellInstance>& cell_instances() const { return cell_instances_; }

  const std::unordered_set<int32_t>& cells() const { return cells_; }

  void set_cell_instances(gsl::span<CellInstance> instances);

private:
  //----------------------------------------------------------------------------
  // Data members

  //! The indices of the cells binned by this filter.
  vector<CellInstance> cell_instances_;

  //! The set of cells used in this filter
  std::unordered_set<int32_t> cells_;

  //! A map from cell/instance indices to filter bin indices.
  std::unordered_map<CellInstance, gsl::index, CellInstanceHash> map_;

  //! Indicates if filter uses only material-filled cells
  bool material_cells_only_;
};

} // namespace openmc

#endif // OPENMC_TALLIES_FILTER_CELL_INSTANCE_H
