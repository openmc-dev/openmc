#include "openmc/tallies/filter_energy.h"

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
  bins_ = get_node_array<double>(node, "bins");
  n_bins_ = bins_.size() - 1;

  // In MG mode, check if the filter bins match the transport bins.
  // We can save tallying time if we know that the tally bins match the energy
  // group structure.  In that case, the matching bin index is simply the group
  // (after flipping for the different ordering of the library and tallying
  // systems).
  if (!settings::run_CE) {
    if (n_bins_ == data::num_energy_groups) {
      matches_transport_groups_ = true;
      for (auto i = 0; i < n_bins_ + 1; i++) {
        if (data::rev_energy_bins[i] != bins_[i]) {
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
      match.bins_.push_back(data::num_energy_groups - p->g + 1);
    } else {
      //TODO: off-by-one
      match.bins_.push_back(data::num_energy_groups - p->last_g + 1);
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
  //TODO: off-by-one
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
    match.bins_.push_back(data::num_energy_groups - p->g + 1);
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
  //TODO: off-by-one
  out << "Outgoing Energy [" << bins_[bin-1] << ", " << bins_[bin] << ")";
  return out.str();
}

//==============================================================================
// Fortran interoperability
//==============================================================================

extern "C" bool energy_filter_matches_transport_groups(EnergyFilter* filt)
{return filt->matches_transport_groups_;}

extern "C" int
energy_filter_search(const EnergyFilter* filt, double val)
{
  if (val < filt->bins_.front() || val > filt->bins_.back()) {
    return -1;
  } else {
    //TODO: off-by-one
    return lower_bound_index(filt->bins_.begin(), filt->bins_.end(), val) + 1;
  }
}

//==============================================================================
// C-API functions
//==============================================================================

extern"C" int
openmc_energy_filter_get_bins(int32_t index, double** energies, int32_t* n)
{
  // Make sure this is a valid index to an allocated filter.
  int err = verify_filter(index);
  if (err) return err;

  // Get a pointer to the filter and downcast.
  auto* filt_base = filter_from_f(index);
  auto* filt = dynamic_cast<EnergyFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to get energy bins on a non-energy filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Output the bins.
  *energies = filt->bins_.data();
  *n = filt->bins_.size();
  return 0;
}

extern "C" int
openmc_energy_filter_set_bins(int32_t index, int32_t n, const double* energies)
{
  // Make sure this is a valid index to an allocated filter.
  int err = verify_filter(index);
  if (err) return err;

  // Get a pointer to the filter and downcast.
  auto* filt_base = filter_from_f(index);
  auto* filt = dynamic_cast<EnergyFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to set energy bins on a non-energy filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Update the filter.
  filt->bins_.clear();
  filt->bins_.resize(n);
  for (int i = 0; i < n; i++) filt->bins_[i] = energies[i];
  filt->n_bins_ = n - 1;
  filter_update_n_bins(index);
  return 0;
}

}// namespace openmc
