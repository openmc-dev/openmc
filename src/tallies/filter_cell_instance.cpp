#include "openmc/tallies/filter_cell_instance.h"

#include <string>

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/xml_interface.h"

namespace openmc {

CellInstanceFilter::CellInstanceFilter(gsl::span<CellInstance> instances)
{
  this->set_cell_instances(instances);
}

void
CellInstanceFilter::from_xml(pugi::xml_node node)
{
  // Get cell IDs/instances
  auto cells = get_node_array<int32_t>(node, "bins");
  Expects(cells.size() % 2 == 0);

  // Convert into vector of CellInstance
  std::vector<CellInstance> instances;
  for (gsl::index i = 0; i < cells.size() / 2; ++i) {
    int32_t cell_id = cells[2*i];
    gsl::index instance = cells[2*i + 1];
    auto search = model::cell_map.find(cell_id);
    if (search == model::cell_map.end()) {
      throw std::runtime_error{fmt::format(
        "Could not find cell {} specified on tally filter.", cell_id)};
    }
    gsl::index index = search->second;
    instances.push_back({index, instance});
  }

  this->set_cell_instances(instances);
}

void
CellInstanceFilter::set_cell_instances(gsl::span<CellInstance> instances)
{
  // Clear existing cells
  cell_instances_.clear();
  cell_instances_.reserve(instances.size());
  map_.clear();

  // Update cells and mapping
  for (auto& x : instances) {
    Expects(x.index_cell >= 0);
    Expects(x.index_cell < model::cells.size());
    const auto& c {model::cells[x.index_cell]};
    if (c->type_ != Fill::MATERIAL) {
      throw std::invalid_argument{fmt::format(
        "Cell {} is not filled with a material. Only material cells can be "
        "used in a cell instance filter.", c->id_)};
    }
    cell_instances_.push_back(x);
    map_[x] = cell_instances_.size() - 1;
  }

  n_bins_ = cell_instances_.size();
}

void
CellInstanceFilter::get_all_bins(const Particle* p, TallyEstimator estimator,
                         FilterMatch& match) const
{
  gsl::index index_cell = p->coord_[p->n_coord_ - 1].cell;
  gsl::index instance = p->cell_instance_;
  auto search = map_.find({index_cell, instance});
  if (search != map_.end()) {
    int index_bin = search->second;
    match.bins_.push_back(index_bin);
    match.weights_.push_back(1.0);
  }
}

void
CellInstanceFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  size_t n = cell_instances_.size();
  xt::xtensor<size_t, 2> data({n, 2});
  for (gsl::index i = 0; i < n; ++i) {
    const auto& x = cell_instances_[i];
    data(i, 0) = model::cells[x.index_cell]->id_;
    data(i, 1) = x.instance;
  }
  write_dataset(filter_group, "bins", data);
}

std::string
CellInstanceFilter::text_label(int bin) const
{
  const auto& x = cell_instances_[bin];
  auto cell_id = model::cells[x.index_cell]->id_;
  return "Cell " + std::to_string(cell_id) + ", Instance "
    + std::to_string(x.instance);
}

} // namespace openmc
