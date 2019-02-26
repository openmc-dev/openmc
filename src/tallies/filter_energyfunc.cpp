#include "openmc/tallies/filter_energyfunc.h"

#include <iomanip>  // for setprecision
#include <ios>  // for scientific
#include <sstream>

#include "openmc/error.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
EnergyFunctionFilter::from_xml(pugi::xml_node node)
{
  if (!settings::run_CE)
    fatal_error("EnergyFunction filters are only supported for "
                "continuous-energy transport calculations");

  if (!check_for_node(node, "energy"))
    fatal_error("Energy grid not specified for EnergyFunction filter.");

  energy_ = get_node_array<double>(node, "energy");

  if (!check_for_node(node, "y"))
    fatal_error("y values not specified for EnergyFunction filter.");

  y_ = get_node_array<double>(node, "y");
}

void
EnergyFunctionFilter::get_all_bins(const Particle* p, int estimator,
                                   FilterMatch& match) const
{
  if (p->last_E >= energy_.front() && p->last_E <= energy_.back()) {
    // Search for the incoming energy bin.
    auto i = lower_bound_index(energy_.begin(), energy_.end(), p->last_E);

    // Compute the interpolation factor between the nearest bins.
    double f = (p->last_E - energy_[i]) / (energy_[i+1] - energy_[i]);

    // Interpolate on the lin-lin grid.
    match.bins_.push_back(0);
    match.weights_.push_back((1-f) * y_[i] + f * y_[i+1]);
  }
}

void
EnergyFunctionFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "energy", energy_);
  write_dataset(filter_group, "y", y_);
}

std::string
EnergyFunctionFilter::text_label(int bin) const
{
  std::stringstream out;
  out << std::scientific << std::setprecision(1)
      << "Energy Function f"
      << "([ " << energy_.front() << ", ..., " << energy_.back() << "]) = "
      << "[" << y_.front() << ", ..., " << y_.back() << "]";
  return out.str();
}

} // namespace openmc
