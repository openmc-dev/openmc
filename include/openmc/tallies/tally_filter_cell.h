#ifndef OPENMC_TALLY_FILTER_CELL_H
#define OPENMC_TALLY_FILTER_CELL_H

#include <cstdint>
#include <sstream>
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

    //TODO: mapping
  }

protected:
  std::vector<int32_t> cells_;
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_CELL_H
