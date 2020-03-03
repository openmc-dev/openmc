#include "openmc/tallies/filter_cell.h"

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
CellFilter::from_xml(pugi::xml_node node)
{
  // Get cell IDs and convert to indices into the global cells vector
  auto cells = get_node_array<int32_t>(node, "bins");
  for (auto& c : cells) {
    auto search = model::cell_map.find(c);
    if (search == model::cell_map.end()) {
      throw std::runtime_error{fmt::format(
        "Could not find cell {} specified on tally filter.", c)};
    }
    c = search->second;
  }

  this->set_cells(cells);
}

void
CellFilter::set_cells(gsl::span<int32_t> cells)
{
  // Clear existing cells
  cells_.clear();
  cells_.reserve(cells.size());
  map_.clear();

  // Update cells and mapping
  for (auto& index : cells) {
    Expects(index >= 0);
    Expects(index < model::cells.size());
    cells_.push_back(index);
    map_[index] = cells_.size() - 1;
  }

  n_bins_ = cells_.size();
}

void
CellFilter::get_all_bins(const Particle* p, TallyEstimator estimator,
                         FilterMatch& match) const
{
  for (int i = 0; i < p->n_coord_; i++) {
    auto search = map_.find(p->coord_[i].cell);
    if (search != map_.end()) {
      match.bins_.push_back(search->second);
      match.weights_.push_back(1.0);
    }
  }
}

void
CellFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  std::vector<int32_t> cell_ids;
  for (auto c : cells_) cell_ids.push_back(model::cells[c]->id_);
  write_dataset(filter_group, "bins", cell_ids);
}

std::string
CellFilter::text_label(int bin) const
{
  return fmt::format("Cell {}", model::cells[cells_[bin]]->id_);
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_cell_filter_get_bins(int32_t index, const int32_t** cells, int32_t* n)
{
  if (int err = verify_filter(index)) return err;

  const auto& filt = model::tally_filters[index].get();
  if (filt->type() != "cell") {
    set_errmsg("Tried to get cells from a non-cell filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  auto cell_filt = static_cast<CellFilter*>(filt);
  *cells = cell_filt->cells().data();
  *n = cell_filt->cells().size();
  return 0;
}

} // namespace openmc
