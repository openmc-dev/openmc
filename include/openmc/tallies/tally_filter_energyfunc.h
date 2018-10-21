#ifndef OPENMC_TALLY_FILTER_ENERGYFUNC_H
#define OPENMC_TALLY_FILTER_ENERGYFUNC_H

#include <iomanip>  // for setprecision
#include <ios>  // for scientific
#include <sstream>
#include <vector>

#include "openmc/error.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/tallies/tally_filter.h"


namespace openmc {

//==============================================================================
//! Multiplies tally scores by an arbitrary function of incident energy
//! described by a piecewise linear-linear interpolation.
//==============================================================================

class EnergyFunctionFilter : public TallyFilter
{
public:
  std::string type() const override {return "energyfunction";}

  EnergyFunctionFilter()
    : TallyFilter {}
  {
    n_bins_ = 1;
  }

  ~EnergyFunctionFilter() = default;

  void
  from_xml(pugi::xml_node node) override
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
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
  {
    if (p->last_E >= energy_.front() && p->last_E <= energy_.back()) {
      // Search for the incoming energy bin.
      auto i = lower_bound_index(energy_.begin(), energy_.end(), p->last_E);

      // Compute the interpolation factor between the nearest bins.
      double f = (p->last_E - energy_[i]) / (energy_[i+1] - energy_[i]);

      // Interpolate on the lin-lin grid.
      match.bins_.push_back(1);
      match.weights_.push_back((1-f) * y_[i] + f * y_[i+1]);
    }
  }

  void
  to_statepoint(hid_t filter_group) const override
  {
    TallyFilter::to_statepoint(filter_group);
    write_dataset(filter_group, "energy", energy_);
    write_dataset(filter_group, "y", y_);
  }

  std::string
  text_label(int bin) const override
  {
    std::stringstream out;
    out << std::scientific << std::setprecision(1)
        << "Energy Function f"
        << "([ " << energy_.front() << ", ..., " << energy_.back() << "]) = "
        << "[" << y_.front() << ", ..., " << y_.back() << "]";
    return out.str();
  }

  std::vector<double> energy_;
  std::vector<double> y_;
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_ENERGYFUNC_H
