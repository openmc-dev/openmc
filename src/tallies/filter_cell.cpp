#include "openmc/tallies/filter_cell.h"

#include <sstream>

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
CellFilter::from_xml(pugi::xml_node node)
{
  cells_ = get_node_array<int32_t>(node, "bins");
  n_bins_ = cells_.size();
}

void
CellFilter::initialize()
{
  // Convert cell IDs to indices of the global array.
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

  // Populate the index->bin map.
  for (int i = 0; i < cells_.size(); i++) {
    map_[cells_[i]] = i;
  }
}

void
CellFilter::get_all_bins(const Particle* p, int estimator,
                         FilterMatch& match) const
{
  for (int i = 0; i < p->n_coord; i++) {
    auto search = map_.find(p->coord[i].cell);
    if (search != map_.end()) {
      //TODO: off-by-one
      match.bins_.push_back(search->second + 1);
      match.weights_.push_back(1.0);
    }
  }
}

void
CellFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  std::vector<int32_t> cell_ids;
  for (auto c : cells_) cell_ids.push_back(cells[c]->id_);
  write_dataset(filter_group, "bins", cell_ids);
}

std::string
CellFilter::text_label(int bin) const
{
  //TODO: off-by-one
  return "Cell " + std::to_string(cells[cells_[bin-1]]->id_);
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_cell_filter_get_bins(int32_t index, int32_t** cells, int32_t* n)
{
  int err = verify_filter(index);
  if (err) return err;

  auto filt = filter_from_f(index);
  if (filt->type() != "cell") {
    set_errmsg("Tried to get cells from a non-cell filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  auto cell_filt = static_cast<CellFilter*>(filt);
  *cells = cell_filt->cells_.data();
  *n = cell_filt->cells_.size();
  return 0;
}

} // namespace openmc
