#include "openmc/tallies/filter_energyfunc.h"

#include <fmt/core.h>

#include "openmc/error.h"
#include "openmc/interpolate.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/xml_interface.h"

namespace openmc {

void EnergyFunctionFilter::from_xml(pugi::xml_node node)
{
  if (!settings::run_CE)
    fatal_error("EnergyFunction filters are only supported for "
                "continuous-energy transport calculations");

  if (!check_for_node(node, "energy"))
    fatal_error("Energy grid not specified for EnergyFunction filter.");

  auto energy = get_node_array<double>(node, "energy");

  if (!check_for_node(node, "y"))
    fatal_error("y values not specified for EnergyFunction filter.");

  auto y = get_node_array<double>(node, "y");
  this->set_data(energy, y);

  // default to linear-linear interpolation
  interpolation_ = Interpolation::lin_lin;
  if (check_for_node(node, "interpolation")) {
    std::string interpolation = get_node_value(node, "interpolation");
    this->set_interpolation(interpolation);
  }
}

void EnergyFunctionFilter::set_data(
  span<const double> energy, span<const double> y)
{
  // Check for consistent sizes with new data
  if (energy.size() != y.size()) {
    fatal_error("Energy grid and y values are not consistent");
  }
  energy_.clear();
  energy_.reserve(energy.size());
  y_.clear();
  y_.reserve(y.size());

  // Copy over energy values, ensuring they are valid
  for (int64_t i = 0; i < energy.size(); ++i) {
    if (i > 0 && energy[i] <= energy[i - 1]) {
      throw std::runtime_error {
        "Energy bins must be monotonically increasing."};
    }
    energy_.push_back(energy[i]);
    y_.push_back(y[i]);
  }
}

void EnergyFunctionFilter::set_interpolation(const std::string& interpolation)
{
  if (interpolation == "histogram") {
    interpolation_ = Interpolation::histogram;
  } else if (interpolation == "linear-linear") {
    interpolation_ = Interpolation::lin_lin;
  } else if (interpolation == "linear-log") {
    interpolation_ = Interpolation::lin_log;
  } else if (interpolation == "log-linear") {
    interpolation_ = Interpolation::log_lin;
  } else if (interpolation == "log-log") {
    interpolation_ = Interpolation::log_log;
  } else if (interpolation == "quadratic") {
    if (energy_.size() < 3)
      fatal_error(
        fmt::format("Quadratic interpolation on EnergyFunctionFilter {} "
                    "requires at least 3 data points.",
          this->id()));
    interpolation_ = Interpolation::quadratic;
  } else if (interpolation == "cubic") {
    if (energy_.size() < 4)
      fatal_error(fmt::format("Cubic interpolation on EnergyFunctionFilter "
                              "{} requires at least 4 data points.",
        this->id()));
    interpolation_ = Interpolation::cubic;
  } else {
    fatal_error(fmt::format(
      "Found invalid interpolation type '{}' on EnergyFunctionFilter {}.",
      interpolation, this->id()));
  }
}

void EnergyFunctionFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  if (p.E_last() >= energy_.front() && p.E_last() <= energy_.back()) {

    double w = interpolate(energy_, y_, p.E_last(), interpolation_);

    // Interpolate on the lin-lin grid.
    match.bins_.push_back(0);
    match.weights_.push_back(w);
  }
}

void EnergyFunctionFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "energy", energy_);
  write_dataset(filter_group, "y", y_);
  hid_t y_dataset = open_dataset(filter_group, "y");
  write_attribute<int>(
    y_dataset, "interpolation", static_cast<int>(interpolation_));
  close_dataset(y_dataset);
}

std::string EnergyFunctionFilter::text_label(int bin) const
{
  return fmt::format(
    "Energy Function f([{:.1e}, ..., {:.1e}]) = [{:.1e}, ..., {:.1e}]",
    energy_.front(), energy_.back(), y_.front(), y_.back());
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int openmc_energyfunc_filter_set_data(
  int32_t index, size_t n, const double* energy, const double* y)
{
  // Ensure this is a valid index to allocated filter
  if (int err = verify_filter(index))
    return err;

  // Get a pointer to the filter
  const auto& filt_base = model::tally_filters[index].get();
  // Downcast to EnergyFunctionFilter
  auto* filt = dynamic_cast<EnergyFunctionFilter*>(filt_base);

  // Check if a valid filter was produced
  if (!filt) {
    set_errmsg(
      "Tried to set interpolation data for non-energy function filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  filt->set_data({energy, n}, {y, n});
  return 0;
}

extern "C" int openmc_energyfunc_filter_get_energy(
  int32_t index, size_t* n, const double** energy)
{
  // ensure this is a valid index to allocated filter
  if (int err = verify_filter(index))
    return err;

  // get a pointer to the filter
  const auto& filt_base = model::tally_filters[index].get();
  // downcast to EnergyFunctionFilter
  auto* filt = dynamic_cast<EnergyFunctionFilter*>(filt_base);

  // check if a valid filter was produced
  if (!filt) {
    set_errmsg(
      "Tried to set interpolation data for non-energy function filter.");
    return OPENMC_E_INVALID_TYPE;
  }
  *energy = filt->energy().data();
  *n = filt->energy().size();
  return 0;
}

extern "C" int openmc_energyfunc_filter_get_y(
  int32_t index, size_t* n, const double** y)
{
  // ensure this is a valid index to allocated filter
  if (int err = verify_filter(index))
    return err;

  // get a pointer to the filter
  const auto& filt_base = model::tally_filters[index].get();
  // downcast to EnergyFunctionFilter
  auto* filt = dynamic_cast<EnergyFunctionFilter*>(filt_base);

  // check if a valid filter was produced
  if (!filt) {
    set_errmsg(
      "Tried to set interpolation data for non-energy function filter.");
    return OPENMC_E_INVALID_TYPE;
  }
  *y = filt->y().data();
  *n = filt->y().size();
  return 0;
}

extern "C" int openmc_energyfunc_filter_set_interpolation(
  int32_t index, const char* interp)
{
  // ensure this is a valid index to allocated filter
  if (int err = verify_filter(index))
    return err;

  // get a pointer to the filter
  const auto& filt_base = model::tally_filters[index].get();
  // downcast to EnergyFunctionFilter
  auto* filt = dynamic_cast<EnergyFunctionFilter*>(filt_base);

  // check if a valid filter was produced
  if (!filt) {
    set_errmsg(
      "Tried to set interpolation data for non-energy function filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Set interpolation
  filt->set_interpolation(interp);
  return 0;
}

extern "C" int openmc_energyfunc_filter_get_interpolation(
  int32_t index, int* interp)
{
  // ensure this is a valid index to allocated filter
  if (int err = verify_filter(index))
    return err;

  // get a pointer to the filter
  const auto& filt_base = model::tally_filters[index].get();
  // downcast to EnergyFunctionFilter
  auto* filt = dynamic_cast<EnergyFunctionFilter*>(filt_base);

  // check if a valid filter was produced
  if (!filt) {
    set_errmsg(
      "Tried to set interpolation data for non-energy function filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  *interp = static_cast<int>(filt->interpolation());
  return 0;
}

} // namespace openmc
