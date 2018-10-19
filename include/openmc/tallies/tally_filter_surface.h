#ifndef OPENMC_TALLY_FILTER_SURFACE_H
#define OPENMC_TALLY_FILTER_SURFACE_H

#include <cstdint>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "openmc/error.h"
#include "openmc/surface.h"
#include "openmc/tallies/tally_filter.h"


namespace openmc {

//==============================================================================
//! Specifies which surface particles are crossing
//==============================================================================

class SurfaceFilter : public TallyFilter
{
public:
  virtual std::string type() const override {return "surface";}

  virtual ~SurfaceFilter() override = default;

  virtual void
  from_xml(pugi::xml_node node) override
  {
    surfaces_ = get_node_array<int32_t>(node, "bins");
    n_bins_ = surfaces_.size();
  }

  virtual void
  initialize() override
  {
    for (auto& s : surfaces_) {
      auto search = surface_map.find(s);
      if (search != surface_map.end()) {
        s = search->second;
      } else {
        std::stringstream err_msg;
        err_msg << "Could not find surface " << s
                << " specified on tally filter.";
        fatal_error(err_msg);
      }
    }

    for (int i = 0; i < surfaces_.size(); i++) {
      map_[surfaces_[i]] = i;
    }
  }

  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
  {
    auto search = map_.find(std::abs(p->surface)-1);
    if (search != map_.end()) {
      match.bins_.push_back(search->second + 1);
      if (p->surface < 0) {
        match.weights_.push_back(-1);
      } else {
        match.weights_.push_back(1);
      }
    }
  }

  virtual void
  to_statepoint(hid_t filter_group) const override
  {
    TallyFilter::to_statepoint(filter_group);
    std::vector<int32_t> surface_ids;
    for (auto c : surfaces_) surface_ids.push_back(surfaces[c]->id_);
    write_dataset(filter_group, "bins", surface_ids);
  }

  virtual std::string
  text_label(int bin) const override
  {
    return "Surface " + std::to_string(surfaces[surfaces_[bin-1]]->id_);
  }

  std::vector<int32_t> surfaces_;
  std::unordered_map<int32_t, int> map_;
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_SURFACE_H
