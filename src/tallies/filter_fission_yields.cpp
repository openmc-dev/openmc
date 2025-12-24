#include "openmc/tallies/filter_fission_yields.h"

#include <pugixml.hpp>

#include "openmc/chain.h"
#include "openmc/endf.h"
#include "openmc/error.h"
#include "openmc/nuclide.h"
#include "openmc/xml_interface.h"

namespace openmc {

void FissionYieldsFilter::from_xml(pugi::xml_node node)
{
  if (!settings::run_CE)
    fatal_error("FissionYieldsFilter filters are only supported for "
                "continuous-energy transport calculations");

  if (!check_for_node(node, "bins"))
    fatal_error("Bins not specified for FissionYieldsFilter.");

  bins_ = get_node_array<std::string>(node, "bins");
  n_bins_ = bins_.size();
}

void FissionYieldsFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  if (p.neutron_xs(p.event_nuclide()).fission > 0) {
    auto nuc = data::nuclides[p.event_nuclide()]->name_;
    if (data::chain_nuclide_map.find(nuc) != data::chain_nuclide_map.end()) {
      auto fy = data::chain_nuclides[data::chain_nuclide_map[nuc]]
                  ->fission_yields()
                  ->yields_;
      for (int i = 0; i < bins_.size(); ++i) {
        if (fy.find(bins_[i]) != fy.end()) {
          match.bins_.push_back(i);
          match.weights_.push_back(fy[bins_[i]](p.E_last()));
        }
      }
    }
  }
}

void FissionYieldsFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", bins_);
}

std::string FissionYieldsFilter::text_label(int bin) const
{
  return fmt::format("Fission Yield [{}]", bins_[bin]);
}

} // namespace openmc
