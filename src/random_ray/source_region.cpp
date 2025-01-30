#include "openmc/random_ray/source_region.h"

namespace openmc {

//==============================================================================
// SourceRegion implementation
//==============================================================================
SourceRegion::SourceRegion(int negroups, bool is_linear)
{
  if (settings::run_mode == RunMode::EIGENVALUE) {
    // If in eigenvalue mode, set starting flux to guess of 1
    scalar_flux_old_.assign(negroups, 1.0);
  } else {
    // If in fixed source mode, set starting flux to guess of zero
    // and initialize external source arrays
    scalar_flux_old_.assign(negroups, 0.0);
    external_source_.assign(negroups, 0.0);
  }

  scalar_flux_new_.assign(negroups, 0.0);
  source_.resize(negroups);
  scalar_flux_final_.assign(negroups, 0.0);

  tally_task_.resize(negroups);

  if (is_linear) {
    source_gradients_.resize(negroups);
    flux_moments_old_.resize(negroups);
    flux_moments_new_.resize(negroups);
    flux_moments_t_.resize(negroups);
  }
}

//==============================================================================
// SourceRegionContainer implementation
//==============================================================================
void SourceRegionContainer::push_back(const SourceRegion& sr)
  {
    num_regions_++;

    // Scalar fields
    material_.push_back(sr.material_);
    lock_.push_back(sr.lock_);
    volume_.push_back(sr.volume_);
    volume_t_.push_back(sr.volume_t_);
    volume_naive_.push_back(sr.volume_naive_);
    position_recorded_.push_back(sr.position_recorded_);
    external_source_present_.push_back(sr.external_source_present_);
    position_.push_back(sr.position_);

    // Only store these fields if is_linear_ is true
    if (is_linear_) {
      centroid_.push_back(sr.centroid_);
      centroid_iteration_.push_back(sr.centroid_iteration_);
      centroid_t_.push_back(sr.centroid_t_);
      mom_matrix_.push_back(sr.mom_matrix_);
      mom_matrix_t_.push_back(sr.mom_matrix_t_);
    }

    // Energy-dependent fields
    for (int g = 0; g < negroups_; ++g) {
      scalar_flux_old_.push_back(sr.scalar_flux_old_[g]);
      scalar_flux_new_.push_back(sr.scalar_flux_new_[g]);
      scalar_flux_final_.push_back(sr.scalar_flux_final_[g]);
      source_.push_back(sr.source_[g]);
      external_source_.push_back(sr.external_source_[g]);

      // Only store these fields if is_linear_ is true
      if (is_linear_) {
        source_gradients_.push_back(sr.source_gradients_[g]);
        flux_moments_old_.push_back(sr.flux_moments_old_[g]);
        flux_moments_new_.push_back(sr.flux_moments_new_[g]);
        flux_moments_t_.push_back(sr.flux_moments_t_[g]);
      }

      // Tally tasks
      tally_task_.emplace_back(sr.tally_task_[g]);
    }
  }

  void SourceRegionContainer::assign(int n_source_regions, const SourceRegion& source_region)
  {
    // Clear existing data
    num_regions_ = 0;
    material_.clear();
    lock_.clear();
    volume_.clear();
    volume_t_.clear();
    volume_naive_.clear();
    position_recorded_.clear();
    external_source_present_.clear();
    position_.clear();

    if (is_linear_) {
      centroid_.clear();
      centroid_iteration_.clear();
      centroid_t_.clear();
      mom_matrix_.clear();
      mom_matrix_t_.clear();
    }

    scalar_flux_old_.clear();
    scalar_flux_new_.clear();
    scalar_flux_final_.clear();
    source_.clear();
    external_source_.clear();

    if (is_linear_) {
      source_gradients_.clear();
      flux_moments_old_.clear();
      flux_moments_new_.clear();
      flux_moments_t_.clear();
    }

    tally_task_.clear();

    // Fill with copies of source_region
    for (int i = 0; i < n_source_regions; ++i) {
      push_back(source_region);
    }
  }

} // namespace openmc