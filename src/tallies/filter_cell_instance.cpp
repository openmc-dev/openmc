#include "openmc/tallies/filter_cell_instance.h"

#include <cassert>
#include <string>

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/xml_interface.h"

namespace openmc {

CellInstanceFilter::CellInstanceFilter(span<CellInstance> instances)
{
  this->set_cell_instances(instances);
}

void CellInstanceFilter::from_xml(pugi::xml_node node)
{
  // Get cell IDs/instances
  auto cells = get_node_array<int32_t>(node, "bins");
  assert(cells.size() % 2 == 0);

  // Convert into vector of CellInstance
  vector<CellInstance> instances;
  for (int64_t i = 0; i < cells.size() / 2; ++i) {
    int32_t cell_id = cells[2 * i];
    int64_t instance = cells[2 * i + 1];
    auto search = model::cell_map.find(cell_id);
    if (search == model::cell_map.end()) {
      throw std::runtime_error {fmt::format(
        "Could not find cell {} specified on tally filter.", cell_id)};
    }
    int64_t index = search->second;
    instances.push_back({index, instance});
  }

  this->set_cell_instances(instances);
}

void CellInstanceFilter::set_cell_instances(span<CellInstance> instances)
{
  // Clear existing cells
  cell_instances_.clear();
  cell_instances_.reserve(instances.size());
  cells_.clear();
  map_.clear();

  // Update cells and mapping
  for (auto& x : instances) {
    assert(x.index_cell >= 0);
    assert(x.index_cell < model::cells.size());
    cell_instances_.push_back(x);
    cells_.insert(x.index_cell);
    map_[x] = cell_instances_.size() - 1;
  }

  n_bins_ = cell_instances_.size();

  material_cells_only_ = true;
  for (const auto& cell_inst : cell_instances_) {
    const auto& c = *model::cells[cell_inst.index_cell];
    if (c.type_ == Fill::MATERIAL)
      continue;
    material_cells_only_ = false;
    break;
  }
}

void CellInstanceFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  int64_t index_cell = p.lowest_coord().cell();
  int64_t instance = p.cell_instance();

  if (cells_.count(index_cell) > 0) {
    auto search = map_.find({index_cell, instance});
    if (search != map_.end()) {
      int index_bin = search->second;
      match.bins_.push_back(index_bin);
      match.weights_.push_back(1.0);
    }
  }

  if (material_cells_only_)
    return;

  for (int i = 0; i < p.n_coord() - 1; i++) {
    int64_t index_cell = p.coord(i).cell();
    // if this cell isn't used on the filter, move on
    if (cells_.count(index_cell) == 0)
      continue;

    // if this cell is used in the filter, check the instance as well
    int64_t instance = cell_instance_at_level(p, i);
    auto search = map_.find({index_cell, instance});
    if (search != map_.end()) {
      match.bins_.push_back(search->second);
      match.weights_.push_back(1.0);
    }
  }
}

void CellInstanceFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  size_t n = cell_instances_.size();
  xt::xtensor<size_t, 2> data({n, 2});
  for (int64_t i = 0; i < n; ++i) {
    const auto& x = cell_instances_[i];
    data(i, 0) = model::cells[x.index_cell]->id_;
    data(i, 1) = x.instance;
  }
  write_dataset(filter_group, "bins", data);
}

std::string CellInstanceFilter::text_label(int bin) const
{
  const auto& x = cell_instances_[bin];
  auto cell_id = model::cells[x.index_cell]->id_;
  return "Cell " + std::to_string(cell_id) + ", Instance " +
         std::to_string(x.instance);
}

} // namespace openmc
