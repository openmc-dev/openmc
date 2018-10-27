#include "openmc/tallies/filter_energy.h"

#include "openmc/constants.h"  // For F90_NONE
#include "openmc/mgxs_interface.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// EnergyFilter implementation
//==============================================================================

// Used to grab the rev_energy_bins array defined in Fortran.
extern "C" double* rev_energy_bins_ptr();

void
EnergyFilter::from_xml(pugi::xml_node node)
{
  bins_ = get_node_array<double>(node, "bins");
  n_bins_ = bins_.size() - 1;

  // In MG mode, check if the filter bins match the transport bins.
  // We can save tallying time if we know that the tally bins match the energy
  // group structure.  In that case, the matching bin index is simply the group
  // (after flipping for the different ordering of the library and tallying
  // systems).
  if (!settings::run_CE) {
    if (n_bins_ == num_energy_groups) {
      matches_transport_groups_ = true;
      double* rev_energy_bins = rev_energy_bins_ptr();
      for (auto i = 0; i < n_bins_ + 1; i++) {
        if (rev_energy_bins[i] != bins_[i]) {
          matches_transport_groups_ = false;
          break;
        }
      }
    }
  }
}

void
EnergyFilter::get_all_bins(const Particle* p, int estimator, FilterMatch& match)
const
{
  if (p->g != F90_NONE && matches_transport_groups_) {
    if (estimator == ESTIMATOR_TRACKLENGTH) {
      //TODO: off-by-one
      match.bins_.push_back(num_energy_groups - p->g + 1);
    } else {
      //TODO: off-by-one
      match.bins_.push_back(num_energy_groups - p->last_g + 1);
    }
    match.weights_.push_back(1.0);

  } else {
    // Get the pre-collision energy of the particle.
    auto E = p->last_E;

    // Bin the energy.
    if (E >= bins_.front() && E <= bins_.back()) {
      //TODO: off-by-one
      auto bin = lower_bound_index(bins_.begin(), bins_.end(), E) + 1;
      match.bins_.push_back(bin);
      match.weights_.push_back(1.0);
    }
  }
}

void
EnergyFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", bins_);
}

std::string
EnergyFilter::text_label(int bin) const
{
  std::stringstream out;
  out << "Incoming Energy [" << bins_[bin-1] << ", " << bins_[bin] << ")";
  return out.str();
}

//==============================================================================
// EnergyoutFilter implementation
//==============================================================================

void
EnergyoutFilter::get_all_bins(const Particle* p, int estimator,
                              FilterMatch& match) const
{
  if (p->g != F90_NONE && matches_transport_groups_) {
    match.bins_.push_back(num_energy_groups - p->g + 1);
    match.weights_.push_back(1.0);

  } else {
    if (p->E >= bins_.front() && p->E <= bins_.back()) {
      //TODO: off-by-one
      auto bin = lower_bound_index(bins_.begin(), bins_.end(), p->E) + 1;
      match.bins_.push_back(bin);
      match.weights_.push_back(1.0);
    }
  }
}

std::string
EnergyoutFilter::text_label(int bin) const
{
  std::stringstream out;
  out << "Outgoing Energy [" << bins_[bin-1] << ", " << bins_[bin] << ")";
  return out.str();
}

}// namespace openmc
