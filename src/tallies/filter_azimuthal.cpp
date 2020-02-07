#include "openmc/tallies/filter_azimuthal.h"

#include <cmath>

#include <fmt/core.h>

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/search.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
AzimuthalFilter::from_xml(pugi::xml_node node)
{
  auto bins = get_node_array<double>(node, "bins");

  if (bins.size() == 1) {
    // Allow a user to input a lone number which will mean that you subdivide
    // [-pi,pi) evenly with the input being the number of bins

    int n_angle = bins[0];
    if (n_angle <= 1) throw std::runtime_error{
      "Number of bins for azimuthal filter must be greater than 1."};

    double d_angle = 2.0 * PI / n_angle;
    bins.resize(n_angle + 1);
    for (int i = 0; i < n_angle; i++) bins[i] = -PI + i * d_angle;
    bins[n_angle] = PI;
  }

  this->set_bins(bins);
}

void AzimuthalFilter::set_bins(gsl::span<double> bins)
{
  // Clear existing bins
  bins_.clear();
  bins_.reserve(bins.size());

  // Copy bins, ensuring they are valid
  for (gsl::index i = 0; i < bins.size(); ++i) {
    if (i > 0 && bins[i] <= bins[i-1]) {
      throw std::runtime_error{"Azimuthal bins must be monotonically increasing."};
    }
    bins_.push_back(bins[i]);
  }

  n_bins_ = bins_.size() - 1;
}

void
AzimuthalFilter::get_all_bins(const Particle* p, TallyEstimator estimator,
                              FilterMatch& match) const
{
  double phi;
  if (estimator == TallyEstimator::TRACKLENGTH) {
    phi = std::atan2(p->u().y, p->u().x);
  } else {
    phi = std::atan2(p->u_last_.y, p->u_last_.x);
  }

  if (phi >= bins_.front() && phi <= bins_.back()) {
    auto bin = lower_bound_index(bins_.begin(), bins_.end(), phi);
    match.bins_.push_back(bin);
    match.weights_.push_back(1.0);
  }
}

void
AzimuthalFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", bins_);
}

std::string
AzimuthalFilter::text_label(int bin) const
{
  return fmt::format("Azimuthal Angle [{}, {})", bins_[bin], bins_[bin+1]);
}

} // namespace openmc
