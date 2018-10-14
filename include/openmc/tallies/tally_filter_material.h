#ifndef OPENMC_TALLY_FILTER_MATERIAL_H
#define OPENMC_TALLY_FILTER_MATERIAL_H

#include <cstdint>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "openmc/error.h"
#include "openmc/material.h"
#include "openmc/tallies/tally_filter.h"


namespace openmc {

class MaterialFilter : public TallyFilter
{
public:
  virtual std::string type() const override {return "material";}

  virtual ~MaterialFilter() override = default;

  virtual void
  from_xml(pugi::xml_node node) override
  {
    materials_ = get_node_array<int32_t>(node, "bins");
    n_bins_ = materials_.size();
  }

  virtual void
  initialize() override
  {
    for (auto& m : materials_) {
      auto search = material_map.find(m);
      if (search != material_map.end()) {
        m = search->second;
      } else {
        std::stringstream err_msg;
        err_msg << "Could not find material " << m
                << " specified on tally filter.";
        fatal_error(err_msg);
      }
    }

    for (int i = 0; i < materials_.size(); i++) {
      map_[materials_[i]] = i;
    }
  }

  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
  {
    auto search = map_.find(p->material - 1);
    if (search != map_.end()) {
      // TODO: off-by-one
      match.bins.push_back(search->second + 1);
      match.weights.push_back(1);
    }
  }

  virtual void
  to_statepoint(hid_t filter_group) const override
  {
    TallyFilter::to_statepoint(filter_group);
    std::vector<int32_t> material_ids;
    for (auto c : materials_) material_ids.push_back(materials[c]->id_);
    write_dataset(filter_group, "bins", material_ids);
  }

  virtual std::string
  text_label(int bin) const override
  {
    return "Material " + std::to_string(materials[materials_[bin-1]]->id_);
  }

  std::vector<int32_t> materials_;
  std::unordered_map<int32_t, int> map_;
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_MATERIAL_H
