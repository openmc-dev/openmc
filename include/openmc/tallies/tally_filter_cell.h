#ifndef OPENMC_TALLY_FILTER_CELL_H
#define OPENMC_TALLY_FILTER_CELL_H

#include <cstdint>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/tallies/tally_filter.h"


namespace openmc {

//==============================================================================
//! Specifies which geometric cells tally events reside in.
//==============================================================================

class CellFilter : public TallyFilter
{
public:
  std::string type() const override {return "cell";}

  ~CellFilter() = default;

  void
  from_xml(pugi::xml_node node) override
  {
    cells_ = get_node_array<int32_t>(node, "bins");
    n_bins_ = cells_.size();
  }

  void
  initialize() override
  {
    for (auto& c : cells_) {
      auto search = cell_map.find(c);
      if (search != cell_map.end()) {
        c = search->second;
      } else {
        std::stringstream err_msg;
        err_msg << "Could not find cell " << c << " specified on tally filter.";
        fatal_error(err_msg);
      }
    }

    for (int i = 0; i < cells_.size(); i++) {
      map_[cells_[i]] = i;
    }
  }

  void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
  {
    for (int i = 0; i < p->n_coord; i++) {
      auto search = map_.find(p->coord[i].cell);
      if (search != map_.end()) {
        // TODO: off-by-one
        match.bins_.push_back(search->second + 1);
        match.weights_.push_back(1);
      }
    }
  }

  void
  to_statepoint(hid_t filter_group) const override
  {
    TallyFilter::to_statepoint(filter_group);
    std::vector<int32_t> cell_ids;
    for (auto c : cells_) cell_ids.push_back(cells[c]->id_);
    write_dataset(filter_group, "bins", cell_ids);
  }

  std::string
  text_label(int bin) const override
  {
    return "Cell " + std::to_string(cells[cells_[bin-1]]->id_);
  }

  std::vector<int32_t> cells_;
  std::unordered_map<int32_t, int> map_;
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_CELL_H
