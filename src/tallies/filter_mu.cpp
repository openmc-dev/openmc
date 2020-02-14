#include "openmc/tallies/filter_mu.h"

#include <fmt/core.h>

#include "openmc/error.h"
#include "openmc/search.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
MuFilter::from_xml(pugi::xml_node node)
{
  auto bins = get_node_array<double>(node, "bins");

  if (bins.size() == 1) {
    // Allow a user to input a lone number which will mean that you subdivide
    // [-1,1) evenly with the input being the number of bins

    int n_angle = bins[0];
    if (n_angle <= 1) throw std::runtime_error{
        "Number of bins for mu filter must be greater than 1."};

    double d_angle = 2.0 / n_angle;
    bins.resize(n_angle + 1);
    for (int i = 0; i < n_angle; i++) bins[i] = -1 + i * d_angle;
    bins[n_angle] = 1;
  }

  this->set_bins(bins);
}

void
MuFilter::set_bins(gsl::span<double> bins)
{
  // Clear existing bins
  bins_.clear();
  bins_.reserve(bins.size());

  // Copy bins, ensuring they are valid
  for (gsl::index i = 0; i < bins.size(); ++i) {
    if (i > 0 && bins[i] <= bins[i-1]) {
      throw std::runtime_error{"Mu bins must be monotonically increasing."};
    }
    bins_.push_back(bins[i]);
  }

  n_bins_ = bins_.size() - 1;
}

void
MuFilter::get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
const
{
  if (p->mu_ >= bins_.front() && p->mu_ <= bins_.back()) {
    auto bin = lower_bound_index(bins_.begin(), bins_.end(), p->mu_);
    match.bins_.push_back(bin);
    match.weights_.push_back(1.0);
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
  return fmt::format("Change-in-Angle [{}, {})", bins_[bin], bins_[bin+1]);
}

} // namespace openmc
