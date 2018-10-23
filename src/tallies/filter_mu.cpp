#include "openmc/tallies/filter_mu.h"

#include <sstream>

#include "openmc/error.h"
#include "openmc/search.h"

namespace openmc {

void
MuFilter::from_xml(pugi::xml_node node)
{
  auto bins = get_node_array<double>(node, "bins");

  if (bins.size() > 1) {
    bins_ = bins;

  } else {
    // Allow a user to input a lone number which will mean that you subdivide
    // [-1,1) evenly with the input being the number of bins

    int n_angle = bins[0];

    if (n_angle <= 1) fatal_error("Number of bins for mu filter must "
                                  "be greater than 1.");

    double d_angle = 2.0 / n_angle;
    bins_.resize(n_angle + 1);
    for (int i = 0; i < n_angle; i++) bins_[i] = -1 + i * d_angle;
    bins_[n_angle] = 1;
  }

  n_bins_ = bins_.size() - 1;
}

void
MuFilter::get_all_bins(Particle* p, int estimator, FilterMatch& match) const
{
  if (p->mu >= bins_[0] && p->mu <= bins_.back()) {
    auto bin = lower_bound_index(bins_.begin(), bins_.end(), p->mu) + 1;
    match.bins_.push_back(bin);
    match.weights_.push_back(1);
  }
}

void
MuFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", bins_);
}

std::string
MuFilter::text_label(int bin) const
{
  std::stringstream out;
  out << "Change-in-Angle [" << bins_[bin-1] << ", " << bins_[bin] << ")";
  return out.str();
}

} // namespace openmc
