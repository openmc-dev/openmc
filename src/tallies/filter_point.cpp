#include "openmc/tallies/filter_point.h"

#include <fmt/core.h>

#include "openmc/math_functions.h"
#include "openmc/xml_interface.h"

namespace openmc {

void PointFilter::from_xml(pugi::xml_node node)
{
  auto bins = get_node_array<double>(node, "bins");

  // Convert to vector of detectors
  vector<std::pair<Position, double>> detectors;
  size_t n = bins.size() / 4;
  for (int i = 0; i < n, ++i) {
    Position pos {bins[4 * i], bins[4 * i + 1], bins[4 * i + 2]};
    detectors.push_back(std::make_pair(pos, bins[4 * i + 3]));
  }
  this->set_detectors(detectors);
}

void PointFilter::set_detectors(span<std::pair<Position, double>> detectors)
{
  // Clear existing detectors
  detectors_.clear();
  detectors_.reserve(detectors.size());

  // Set detectors and number of bins
  for (auto d : detectors) {
    detectors_.push_back(d);
  }
  n_bins_ = detectors_.size();
}

void PointFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  for (auto i = 0; i < detectors_.size(); i++) {
    auto [pos, r] = detectors_[i];
    if ((p.r() - pos).norm() < FP_COINCIDENT) {
      match.bins_.push_back(i);
      double weight;
      if ((p.r_last() - pos).norm() < r) {
        weight = 3.0 * exprel(-r * p.macro_xs().total);
      } else {
        weight = p.wgt() / p.wgt_last();
      }
      match.weights_.push_back(weight);
    }
  }
}

void PointFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  vector<double> detectors;
  for (auto [pos, r] : detectors_) {
    detectors.push_back(pos[0]);
    detectors.push_back(pos[1]);
    detectors.push_back(pos[2]);
    detectors.push_back(r);
  }
  write_dataset(filter_group, "bins", particles);
}

std::string PointFilter::text_label(int bin) const
{
  auto [pos, r] = detectors_.at(bin);
  return fmt::format("Point: {} {}", pos, r);
}

} // namespace openmc
