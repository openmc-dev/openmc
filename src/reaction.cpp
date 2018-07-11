#include "reaction.h"

#include <string>
#include <utility> // for move

#include "hdf5_interface.h"

namespace openmc {

Reaction::Reaction(hid_t group, const std::vector<int>& temperatures)
{
  read_attribute(group, "Q_value", q_value_);
  read_attribute(group, "mt", mt_);
  int cm;
  read_attribute(group, "center_of_mass", cm);
  scatter_in_cm_ = (cm == 1);

  // Read cross section and threshold_idx data
  for (auto t : temperatures) {
    // Get group corresponding to temperature
    std::string temp_str {std::to_string(t) + "K"};
    hid_t temp_group = open_group(group, temp_str.c_str());
    hid_t dset = open_dataset(temp_group, "xs");

    // Get threshold index
    TemperatureXS xs;
    read_attribute(dset, "threshold_idx", xs.threshold);

    // Read cross section values
    read_dataset(dset, xs.value);
    close_dataset(dset);
    close_group(temp_group);

    // create new entry in xs vector
    xs_.push_back(std::move(xs));
  }

  // Read products
  for (const auto& name : group_names(group)) {
    if (name.rfind("product_", 0) == 0) {
      hid_t pgroup = open_group(group, name.c_str());
      products_.emplace_back(pgroup);
      close_group(pgroup);
    }
  }
}

}