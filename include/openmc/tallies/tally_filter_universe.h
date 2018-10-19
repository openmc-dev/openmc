#ifndef OPENMC_TALLY_FILTER_UNIVERSE_H
#define OPENMC_TALLY_FILTER_UNIVERSE_H

#include <cstdint>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/tallies/tally_filter.h"


namespace openmc {

//==============================================================================
//! Specifies which geometric universes tally events reside in.
//==============================================================================

class UniverseFilter : public TallyFilter
{
public:
  virtual std::string type() const override {return "universe";}

  virtual ~UniverseFilter() override = default;

  virtual void
  from_xml(pugi::xml_node node) override
  {
    universes_ = get_node_array<int32_t>(node, "bins");
    n_bins_ = universes_.size();
  }

  virtual void
  initialize() override
  {
    for (auto& u : universes_) {
      auto search = universe_map.find(u);
      if (search != universe_map.end()) {
        u = search->second;
      } else {
        std::stringstream err_msg;
        err_msg << "Could not find universe " << u
                << " specified on tally filter.";
        fatal_error(err_msg);
      }
    }

    for (int i = 0; i < universes_.size(); i++) {
      map_[universes_[i]] = i;
    }
  }

  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
  {
    for (int i = 0; i < p->n_coord; i++) {
      auto search = map_.find(p->coord[i].universe);
      if (search != map_.end()) {
        match.bins_.push_back(search->second + 1);
        match.weights_.push_back(1);
      }
    }
  }

  virtual void
  to_statepoint(hid_t filter_group) const override
  {
    TallyFilter::to_statepoint(filter_group);
    std::vector<int32_t> universe_ids;
    for (auto u : universes_) universe_ids.push_back(universes[u]->id_);
    write_dataset(filter_group, "bins", universe_ids);
  }

  virtual std::string
  text_label(int bin) const override
  {
    return "Universe " + std::to_string(universes[universes_[bin-1]]->id_);
  }

  std::vector<int32_t> universes_;
  std::unordered_map<int32_t, int> map_;
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_UNIVERSE_H
