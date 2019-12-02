#ifndef OPENMC_TALLIES_FILTER_CELL_INSTANCE_H
#define OPENMC_TALLIES_FILTER_CELL_INSTANCE_H

#include <cstdint>
#include <functional> // for hash
#include <unordered_map>
#include <vector>

#include <gsl/gsl>

#include "openmc/tallies/filter.h"

//==============================================================================
//! Define an instance of a particular cell
//==============================================================================

namespace openmc {
struct CellInstance {
  //! Check for equality
  bool operator==(const CellInstance& other) const
  { return index_cell == other.index_cell && instance == other.instance; }

  gsl::index index_cell;
  gsl::index instance;
};
} // namespace openmc

namespace std {
template<>
struct hash<openmc::CellInstance> {
  // Taken from https://stackoverflow.com/a/17017281
  std::size_t operator()(const openmc::CellInstance& k) const
  {
    std::size_t res = 17;
    res = 31 * res + std::hash<gsl::index>()(k.index_cell);
    res = 31 * res + std::hash<gsl::index>()(k.instance);
    return res;
  }
};
} // namespace std


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

  std::string type() const override {return "cellinstance";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const std::vector<CellInstance>& cell_instances() const { return cell_instances_; }

  void set_cell_instances(gsl::span<CellInstance> instances);

private:
  //----------------------------------------------------------------------------
  // Data members

  //! The indices of the cells binned by this filter.
  std::vector<CellInstance> cell_instances_;

  //! A map from cell/instance indices to filter bin indices.
  std::unordered_map<CellInstance, gsl::index> map_;
};

} // namespace openmc

#endif // OPENMC_TALLIES_FILTER_CELL_INSTANCE_H
