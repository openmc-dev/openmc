#include "openmc/tallies/filter_collision.h"

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/constants.h"  // For F90_NONE
#include "openmc/mgxs_interface.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// CollisionFilter implementation
//==============================================================================

void
CollisionFilter::from_xml(pugi::xml_node node)
{
  auto bins = get_node_array<int>(node, "bins");
  this->set_bins(bins);
}

void
CollisionFilter::set_bins(gsl::span<const int> bins)
{
  // Clear existing bins
  bins_.clear();
  bins_.reserve(bins.size());

  // Copy bins, ensuring they are valid
  for (gsl::index i = 0; i < bins.size(); ++i) {
    if (i > 0 && bins[i] <= bins[i-1]) {
      throw std::runtime_error{"Number of Collisions bins must be monotonically increasing."};
    }
    bins_.push_back(bins[i]);
  }

  n_bins_ = bins_.size();

}

void
CollisionFilter::get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match)
const
{
    // Get the number of collisions for the particle
    auto n = p.n_collision_;

    // Bin the collision number. Must fit exactly the desired collision number .
    if (n >= bins_.front() && n <= bins_.back()) {
      auto it = find(bins_.begin(), bins_.end(), n);
      if (it != bins_.end()){
        size_t bin = it - bins_.begin();
        if (int(bins_[bin]) == n){ 
          match.bins_.push_back(bin);
          match.weights_.push_back(1.0);
        }
      }
    }
}

void
CollisionFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", bins_);
}

std::string
CollisionFilter::text_label(int bin) const
{
  return fmt::format("Collision Number {}", int(bins_[bin]));
}

//==============================================================================
// C-API functions
//==============================================================================

extern"C" int
openmc_collision_filter_get_bins(int32_t index, const int** energies, size_t* n)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<CollisionFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to get collision bins on a non-collision filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Output the bins.
  *energies = filt->bins().data();
  *n = filt->bins().size();
  return 0;
}

extern "C" int
openmc_collision_filter_set_bins(int32_t index, size_t n, const int* energies)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<CollisionFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to set collision bins on a non-collision filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Update the filter.
  filt->set_bins({energies, n});
  return 0;
}

}// namespace openmc
