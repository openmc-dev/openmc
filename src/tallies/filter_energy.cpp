#include "openmc/tallies/filter_energy.h"

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/constants.h"  // For F90_NONE
#include "openmc/mgxs_interface.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// EnergyFilter implementation
//==============================================================================

void
EnergyFilter::from_xml(pugi::xml_node node)
{
  auto bins = get_node_array<double>(node, "bins");
  this->set_bins(bins);
}

void
EnergyFilter::set_bins(gsl::span<const double> bins)
{
  // Clear existing bins
  bins_.clear();
  bins_.reserve(bins.size());

  // Copy bins, ensuring they are valid
  for (gsl::index i = 0; i < bins.size(); ++i) {
    if (i > 0 && bins[i] <= bins[i-1]) {
      throw std::runtime_error{"Energy bins must be monotonically increasing."};
    }
    bins_.push_back(bins[i]);
  }

  n_bins_ = bins_.size() - 1;

  // In MG mode, check if the filter bins match the transport bins.
  // We can save tallying time if we know that the tally bins match the energy
  // group structure.  In that case, the matching bin index is simply the group
  // (after flipping for the different ordering of the library and tallying
  // systems).
  if (!settings::run_CE) {
    if (n_bins_ == data::mg.num_energy_groups_) {
      matches_transport_groups_ = true;
      for (gsl::index i = 0; i < n_bins_ + 1; ++i) {
        if (data::mg.rev_energy_bins_[i] != bins_[i]) {
          matches_transport_groups_ = false;
          break;
        }
      }
    }
  }
}

void
EnergyFilter::get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
const
{
  if (p->g_ != F90_NONE && matches_transport_groups_) {
    if (estimator == TallyEstimator::TRACKLENGTH) {
      match.bins_.push_back(data::mg.num_energy_groups_ - p->g_ - 1);
    } else {
      match.bins_.push_back(data::mg.num_energy_groups_ - p->g_last_ - 1);
    }
    match.weights_.push_back(1.0);

  } else {
    // Get the pre-collision energy of the particle.
    auto E = p->E_last_;

    // Bin the energy.
    if (E >= bins_.front() && E <= bins_.back()) {
      auto bin = lower_bound_index(bins_.begin(), bins_.end(), E);
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
  return fmt::format("Incoming Energy [{}, {})", bins_[bin], bins_[bin+1]);
}

//==============================================================================
// EnergyoutFilter implementation
//==============================================================================

void
EnergyoutFilter::get_all_bins(const Particle* p, TallyEstimator estimator,
                              FilterMatch& match) const
{
  if (p->g_ != F90_NONE && matches_transport_groups_) {
    match.bins_.push_back(data::mg.num_energy_groups_ - p->g_ - 1);
    match.weights_.push_back(1.0);

  } else {
    if (p->E_ >= bins_.front() && p->E_ <= bins_.back()) {
      auto bin = lower_bound_index(bins_.begin(), bins_.end(), p->E_);
      match.bins_.push_back(bin);
      match.weights_.push_back(1.0);
    }
  }
}

std::string
EnergyoutFilter::text_label(int bin) const
{
  return fmt::format("Outgoing Energy [{}, {})", bins_[bin], bins_[bin+1]);
}

//==============================================================================
// C-API functions
//==============================================================================

extern"C" int
openmc_energy_filter_get_bins(int32_t index, const double** energies, size_t* n)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<EnergyFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to get energy bins on a non-energy filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Output the bins.
  *energies = filt->bins().data();
  *n = filt->bins().size();
  return 0;
}

extern "C" int
openmc_energy_filter_set_bins(int32_t index, size_t n, const double* energies)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<EnergyFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to set energy bins on a non-energy filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Update the filter.
  filt->set_bins({energies, n});
  return 0;
}

}// namespace openmc
