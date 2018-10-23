#include "openmc/tallies/filter_polar.h"

#include <sstream>

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/search.h"

namespace openmc {

void
PolarFilter::from_xml(pugi::xml_node node)
{
  auto bins = get_node_array<double>(node, "bins");

  if (bins.size() > 1) {
    bins_ = bins;

  } else {
    // Allow a user to input a lone number which will mean that you subdivide
    // [0,pi] evenly with the input being the number of bins

    int n_angle = bins[0];

    if (n_angle <= 1) fatal_error("Number of bins for polar filter must "
                                  "be greater than 1.");

    double d_angle = PI / n_angle;
    bins_.resize(n_angle + 1);
    for (int i = 0; i < n_angle; i++) bins_[i] = i * d_angle;
    bins_[n_angle] = PI;
  }

  n_bins_ = bins_.size() - 1;
}

void
PolarFilter::get_all_bins(Particle* p, int estimator, FilterMatch& match)
const
{
  double theta;
  if (estimator == ESTIMATOR_TRACKLENGTH) {
    theta = std::acos(p->coord[0].uvw[2]);
  } else {
    theta = std::acos(p->last_uvw[2]);
  }

  if (theta >= bins_[0] && theta <= bins_.back()) {
    auto bin = lower_bound_index(bins_.begin(), bins_.end(), theta) + 1;
    match.bins_.push_back(bin);
    match.weights_.push_back(1);
  }
}

void
PolarFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", bins_);
}

std::string
PolarFilter::text_label(int bin) const
{
  std::stringstream out;
  out << "Polar Angle [" << bins_[bin-1] << ", " << bins_[bin] << ")";
  return out.str();
}

} // namespace openmc
