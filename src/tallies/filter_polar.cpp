#include "openmc/tallies/filter.h"

#include <fmt/core.h>

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/search.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
Filter::PolarFilter_from_xml(pugi::xml_node node)
{
  auto bins = get_node_array<double>(node, "bins");

  if (bins.size() == 1) {
    // Allow a user to input a lone number which will mean that you subdivide
    // [0,pi] evenly with the input being the number of bins

    int n_angle = bins[0];
    if (n_angle <= 1) throw std::runtime_error{
      "Number of bins for polar filter must be greater than 1."};

    double d_angle = PI / n_angle;
    bins.resize(n_angle + 1);
    for (int i = 0; i < n_angle; i++) bins[i] = i * d_angle;
    bins[n_angle] = PI;
  }

  this->set_bins(bins);
}

void
Filter::PolarFilter_set_bins(gsl::span<const double> bins)
{
  // Clear existing bins
  bins_.clear();
  bins_.reserve(bins.size());

  // Copy bins, ensuring they are valid
  for (gsl::index i = 0; i < bins.size(); ++i) {
    if (i > 0 && bins[i] <= bins[i-1]) {
      throw std::runtime_error{"Polar bins must be monotonically increasing."};
    }
    bins_.push_back(bins[i]);
  }

  n_bins_ = bins_.size() - 1;
}

void
Filter::PolarFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match)
const
{
  double z = (estimator == TallyEstimator::TRACKLENGTH) ? p.u().z : p.u_last_.z;
  double theta = std::acos(z);

  if (theta >= bins_.front() && theta <= bins_.back()) {
    auto bin = lower_bound_index(bins_.begin(), bins_.end(), theta);
    //match.bins_.push_back(bin);
    //match.weights_.push_back(1.0);
    match.push_back(bin, 1.0);
  }
}

void
Filter::PolarFilter_to_statepoint(hid_t filter_group) const
{
  write_dataset(filter_group, "bins", bins_);
}

std::string
Filter::PolarFilter_text_label(int bin) const
{
  return fmt::format("Polar Angle [{}, {})", bins_[bin], bins_[bin+1]);
}

} // namespace openmc
