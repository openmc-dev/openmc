#include "openmc/tallies/filter_timed_mesh.h"

#include <fmt/core.h>
#include <gsl/gsl-lite.hpp>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/mesh.h"
#include "openmc/search.h"
#include "openmc/xml_interface.h"

namespace openmc {

void TimedMeshFilter::from_xml(pugi::xml_node node)
{
  pugi::xml_node bin_node = node.child("bins");

  //----------------------------------------------------------------------------
  // Mesh

  auto mesh_bins_ = get_node_array<int32_t>(bin_node, "mesh_bins");
  if (mesh_bins_.size() != 1) {
    fatal_error(
      "Only one mesh can be specified per " + type_str() + " mesh filter.");
  }

  auto id = mesh_bins_[0];
  auto search = model::mesh_map.find(id);
  if (search != model::mesh_map.end()) {
    set_mesh(search->second);
  } else {
    fatal_error(
      fmt::format("Could not find mesh {} specified on tally filter.", id));
  }

  if (check_for_node(node, "translation")) {
    set_translation(get_node_array<double>(node, "translation"));
  }

  //----------------------------------------------------------------------------
  // Time bins

  auto time_grid = get_node_array<double>(bin_node, "time_bins");
  this->set_time_grid(time_grid);
}

void TimedMeshFilter::set_mesh(int32_t mesh)
{
  // perform any additional perparation for mesh tallies here
  mesh_ = mesh;
  model::meshes[mesh_]->prepare_for_point_location();

  reset_bins();
}

void TimedMeshFilter::set_time_grid(gsl::span<const double> time_grid)
{
  // Clear existing bins
  time_grid_.clear();
  time_grid_.reserve(time_grid.size());

  // Ensure time grid is sorted and don't have duplicates
  if (std::adjacent_find(time_grid.cbegin(), time_grid.cend(),
        std::greater_equal<>()) != time_grid.end()) {
    throw std::runtime_error {"Time grid must be monotonically increasing."};
  }

  // Copy grid
  std::copy(
    time_grid.cbegin(), time_grid.cend(), std::back_inserter(time_grid_));

  reset_bins();
}

void TimedMeshFilter::reset_bins()
{
  mesh_n_bins_ = model::meshes[mesh_]->n_bins();
  n_bins_ = (time_grid_.size() - 1) * mesh_n_bins_;
}

void TimedMeshFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // Get the start/end time of the particle for this track
  const auto t_start = p.time_last();
  const auto t_end = p.time();

  // If time interval is entirely out of time bin range, exit
  if (t_end < time_grid_.front() || t_start >= time_grid_.back())
    return;

  // Get the start/end positions, direction, and speed
  Position last_r = p.r_last();
  Position r = p.r();
  Position u = p.u();
  const auto speed = p.speed();

  // apply translation if present
  if (translated_) {
    last_r -= translation();
    r -= translation();
  }

  if (estimator != TallyEstimator::TRACKLENGTH) {
    // -------------------------------------------------------------------------
    // Non-tracklength estimators
    // Find a match based on the exact time-position of the particle

    // Proceed only if inside time grid
    if (t_end < time_grid_.back())
      return;

    auto time_bin =
      lower_bound_index(time_grid_.begin(), time_grid_.end(), t_end);
    auto mesh_bin = model::meshes[mesh_]->get_bin(r);

    // Inside the mesh?
    if (mesh_bin < 0)
      return;

    auto bin = time_bin * mesh_n_bins_ + mesh_bin;
    match.bins_.push_back(bin);
    match.weights_.push_back(1.0);

  } else {
    // -------------------------------------------------------------------------
    // For tracklength estimator, we have to check the start/end time-position
    // of the current track and find where it overlaps with time-mesh bins and
    // score accordingly.

    // Determine first time bin containing a portion of time interval
    auto i_time_bin =
      lower_bound_index(time_grid_.begin(), time_grid_.end(), t_start);

    // std::cout<<last_r.x<<"  "<<r.x<<"  "<<u.x<<"  "<<t_start<<"  "<<t_end<<"
    // "<<speed<<"\n";

    // If time interval is zero, add a match corresponding to the starting time
    if (t_end == t_start) {
      model::meshes[mesh_]->bins_crossed(
        last_r, r, u, match.bins_, match.weights_);
      // Offset the bin location accordingly
      match.bins_.back() += i_time_bin * mesh_n_bins_;

      /*
      std::cout<<"end = start\n";
      double sum {0.0};
      for (int i = 0; i < match.bins_.size(); i++) {
        std::cout<<"    "<<match.bins_[i]<<"  "<<match.weights_[i]<<"\n";
        sum += match.weights_[i];
      }
      std::cout<<sum<<"\n";
      */

      return;
    }

    // Find matching bins
    double dt_total = t_end - t_start;
    for (; i_time_bin < time_grid_.size() - 1; ++i_time_bin) {
      const double t_left = std::max(t_start, time_grid_[i_time_bin]);
      const double t_right = std::min(t_end, time_grid_[i_time_bin + 1]);

      // Time interval and its fraction
      const double dt = t_right - t_left;
      const double fraction = dt / dt_total;

      // Starting and ending position in this time interval
      const Position r_start = last_r + u * speed * (t_left - t_start);
      const Position r_end = r_start + u * speed * dt;

      // Mesh sweep in this time interval
      const auto n_match_old = match.bins_.size();
      model::meshes[mesh_]->bins_crossed(
        r_start, r_end, u, match.bins_, match.weights_);
      const auto n_match = match.bins_.size() - n_match_old;

      // Update the newly-added bins and weights
      const auto offset = i_time_bin * mesh_n_bins_;
      for (int i = 0; i < n_match; i++) {
        match.bins_[match.bins_.size() - i - 1] += offset;
        match.weights_[match.weights_.size() - i - 1] *= fraction;
      }

      if (t_end < time_grid_[i_time_bin + 1])
        break;
    }

    /*
    double sum {0.0};
    for (int i = 0; i < match.bins_.size(); i++) {
      std::cout<<"    "<<match.bins_[i]<<"  "<<match.weights_[i]<<"\n";
      sum += match.weights_[i];
    }
    std::cout<<sum<<"\n";
    */
  }
}

void TimedMeshFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "time_bins", time_grid_);
  write_dataset(filter_group, "mesh_bins", model::meshes[mesh_]->id_);
  if (translated_) {
    write_dataset(filter_group, "translation", translation_);
  }
}

std::string TimedMeshFilter::text_label(int bin) const
{
  int bin_time = bin / mesh_n_bins_;
  int bin_mesh = bin % mesh_n_bins_;

  auto& mesh = *model::meshes.at(mesh_);

  std::string label_time =
    fmt::format("Time [{}, {}) : ", time_grid_[bin], time_grid_[bin + 1]);
  std::string label_mesh = mesh.bin_label(bin_mesh);

  return label_time + label_mesh;
}

void TimedMeshFilter::set_translation(const Position& translation)
{
  translated_ = true;
  translation_ = translation;
}

void TimedMeshFilter::set_translation(const double translation[3])
{
  this->set_translation({translation[0], translation[1], translation[2]});
}

} // namespace openmc
