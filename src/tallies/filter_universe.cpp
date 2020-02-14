#include "openmc/tallies/filter_universe.h"

#include <fmt/core.h>

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
UniverseFilter::from_xml(pugi::xml_node node)
{
  // Get material IDs and convert to indices in the global materials vector
  auto universes = get_node_array<int32_t>(node, "bins");
  for (auto& u : universes) {
    auto search = model::universe_map.find(u);
    if (search == model::universe_map.end()) {
      throw std::runtime_error{fmt::format(
        "Could not find universe {} specified on tally filter.", u)};
    }
    u = search->second;
  }

  this->set_universes(universes);
}

void
UniverseFilter::set_universes(gsl::span<int32_t> universes)
{
  // Clear existing universes
  universes_.clear();
  universes_.reserve(universes.size());
  map_.clear();

  // Update universes and mapping
  for (auto& index : universes) {
    Expects(index >= 0);
    Expects(index < model::universes.size());
    universes_.push_back(index);
    map_[index] = universes_.size() - 1;
  }

  n_bins_ = universes_.size();
}

void
UniverseFilter::get_all_bins(const Particle* p, TallyEstimator estimator,
                             FilterMatch& match) const
{
  for (int i = 0; i < p->n_coord_; i++) {
    auto search = map_.find(p->coord_[i].universe);
    if (search != map_.end()) {
      match.bins_.push_back(search->second);
      match.weights_.push_back(1.0);
    }
  }
}

void
UniverseFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  std::vector<int32_t> universe_ids;
  for (auto u : universes_) universe_ids.push_back(model::universes[u]->id_);
  write_dataset(filter_group, "bins", universe_ids);
}

std::string
UniverseFilter::text_label(int bin) const
{
  return fmt::format("Universe {}", model::universes[universes_[bin]]->id_);
}

} // namespace openmc
