#include "openmc/tallies/filter_universe.h"

#include <sstream>

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
UniverseFilter::from_xml(pugi::xml_node node)
{
  universes_ = get_node_array<int32_t>(node, "bins");
  n_bins_ = universes_.size();
}

void
UniverseFilter::initialize()
{
  // Convert universe IDs to indices of the global array.
  for (auto& u : universes_) {
    auto search = model::universe_map.find(u);
    if (search != model::universe_map.end()) {
      u = search->second;
    } else {
      std::stringstream err_msg;
      err_msg << "Could not find universe " << u
              << " specified on tally filter.";
      fatal_error(err_msg);
    }
  }

  // Populate the index->bin map.
  for (int i = 0; i < universes_.size(); i++) {
    map_[universes_[i]] = i;
  }
}

void
UniverseFilter::get_all_bins(const Particle* p, int estimator,
                             FilterMatch& match) const
{
  for (int i = 0; i < p->n_coord; i++) {
    auto search = map_.find(p->coord[i].universe);
    if (search != map_.end()) {
      //TODO: off-by-one
      match.bins_.push_back(search->second + 1);
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
  //TODO: off-by-one
  return "Universe " + std::to_string(model::universes[universes_[bin-1]]->id_);
}

} // namespace openmc
