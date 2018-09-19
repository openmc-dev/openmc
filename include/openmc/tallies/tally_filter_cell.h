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

class CellFilter : public TallyFilter
{
public:
  virtual void
  from_xml(pugi::xml_node node)
  {
    cells_ = get_node_array<int32_t>(node, "bins");
  }

  virtual void
  initialize()
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

  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match) const
  {
    for (int i = 0; i < p->n_coord; i++) {
      auto search = map_.find(p->coord[i].cell);
      if (search != map_.end()) {
        // TODO: off-by-one
        match.bins.push_back(search->second + 1);
        match.weights.push_back(1);
      }
    }
  }

protected:
  std::vector<int32_t> cells_;
  std::unordered_map<int32_t, int> map_;
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_CELL_H
