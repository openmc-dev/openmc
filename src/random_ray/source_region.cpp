#include "openmc/random_ray/source_region.h"

#include "openmc/error.h"
#include "openmc/message_passing.h"
#include "openmc/simulation.h"

namespace openmc {

//==============================================================================
// SourceRegionHandle implementation
//==============================================================================
SourceRegionHandle::SourceRegionHandle(SourceRegion& sr)
  : negroups_(sr.scalar_flux_old_.size()), material_(&sr.material_),
    is_small_(&sr.is_small_), n_hits_(&sr.n_hits_),
    is_linear_(sr.source_gradients_.size() > 0), lock_(&sr.lock_),
    volume_(&sr.volume_), volume_t_(&sr.volume_t_), volume_sq_(&sr.volume_sq_),
    volume_sq_t_(&sr.volume_sq_t_), volume_naive_(&sr.volume_naive_),
    position_recorded_(&sr.position_recorded_),
    external_source_present_(&sr.external_source_present_),
    position_(&sr.position_), centroid_(&sr.centroid_),
    centroid_iteration_(&sr.centroid_iteration_), centroid_t_(&sr.centroid_t_),
    mom_matrix_(&sr.mom_matrix_), mom_matrix_t_(&sr.mom_matrix_t_),
    volume_task_(&sr.volume_task_), mesh_(&sr.mesh_),
    parent_sr_(&sr.parent_sr_), scalar_flux_old_(sr.scalar_flux_old_.data()),
    scalar_flux_new_(sr.scalar_flux_new_.data()), source_(sr.source_.data()),
    external_source_(sr.external_source_.data()),
    scalar_flux_final_(sr.scalar_flux_final_.data()),
    source_gradients_(sr.source_gradients_.data()),
    flux_moments_old_(sr.flux_moments_old_.data()),
    flux_moments_new_(sr.flux_moments_new_.data()),
    flux_moments_t_(sr.flux_moments_t_.data()),
    tally_task_(sr.tally_task_.data())
{}

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

SourceRegion::SourceRegion(const SourceRegionHandle& handle, int64_t parent_sr)
  : SourceRegion(handle.negroups_, handle.is_linear_)
{
  material_ = handle.material();
  mesh_ = handle.mesh();
  parent_sr_ = parent_sr;
  for (int g = 0; g < scalar_flux_new_.size(); g++) {
    scalar_flux_old_[g] = handle.scalar_flux_old(g);
    source_[g] = handle.source(g);
  }

  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    external_source_present_ = handle.external_source_present();
    for (int g = 0; g < scalar_flux_new_.size(); g++) {
      external_source_[g] = handle.external_source(g);
    }
  }
}

//==============================================================================
// SourceRegionContainer implementation
//==============================================================================

void SourceRegionContainer::push_back(const SourceRegion& sr)
{
  n_source_regions_++;

  // Scalar fields
  material_.push_back(sr.material_);
  is_small_.push_back(sr.is_small_);
  n_hits_.push_back(sr.n_hits_);
  lock_.push_back(sr.lock_);
  volume_.push_back(sr.volume_);
  volume_t_.push_back(sr.volume_t_);
  volume_sq_.push_back(sr.volume_sq_);
  volume_sq_t_.push_back(sr.volume_sq_t_);
  volume_naive_.push_back(sr.volume_naive_);
  position_recorded_.push_back(sr.position_recorded_);
  external_source_present_.push_back(sr.external_source_present_);
  position_.push_back(sr.position_);
  volume_task_.push_back(sr.volume_task_);
  mesh_.push_back(sr.mesh_);
  parent_sr_.push_back(sr.parent_sr_);

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
    if (settings::run_mode == RunMode::FIXED_SOURCE) {
      external_source_.push_back(sr.external_source_[g]);
    }

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

void SourceRegionContainer::assign(
  int n_source_regions, const SourceRegion& source_region)
{
  // Clear existing data
  n_source_regions_ = 0;
  material_.clear();
  is_small_.clear();
  n_hits_.clear();
  lock_.clear();
  volume_.clear();
  volume_t_.clear();
  volume_sq_.clear();
  volume_sq_t_.clear();
  volume_naive_.clear();
  position_recorded_.clear();
  external_source_present_.clear();
  position_.clear();
  mesh_.clear();
  parent_sr_.clear();

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
  volume_task_.clear();

  // Fill with copies of source_region
  for (int i = 0; i < n_source_regions; ++i) {
    push_back(source_region);
  }
}

void SourceRegionContainer::flux_swap()
{
  scalar_flux_old_.swap(scalar_flux_new_);
  if (is_linear_) {
    flux_moments_old_.swap(flux_moments_new_);
  }
}

SourceRegionHandle SourceRegionContainer::get_source_region_handle(int64_t sr)
{
  SourceRegionHandle handle;
  handle.negroups_ = negroups();
  handle.material_ = &material(sr);
  handle.is_small_ = &is_small(sr);
  handle.n_hits_ = &n_hits(sr);
  handle.is_linear_ = is_linear();
  handle.lock_ = &lock(sr);
  handle.volume_ = &volume(sr);
  handle.volume_t_ = &volume_t(sr);
  handle.volume_sq_ = &volume_sq(sr);
  handle.volume_sq_t_ = &volume_sq_t(sr);
  handle.volume_naive_ = &volume_naive(sr);
  handle.position_recorded_ = &position_recorded(sr);
  handle.external_source_present_ = &external_source_present(sr);
  handle.position_ = &position(sr);
  handle.volume_task_ = &volume_task(sr);
  handle.mesh_ = &mesh(sr);
  handle.parent_sr_ = &parent_sr(sr);
  handle.scalar_flux_old_ = &scalar_flux_old(sr, 0);
  handle.scalar_flux_new_ = &scalar_flux_new(sr, 0);
  handle.source_ = &source(sr, 0);
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    handle.external_source_ = &external_source(sr, 0);
  } else {
    handle.external_source_ = nullptr;
  }
  handle.scalar_flux_final_ = &scalar_flux_final(sr, 0);
  handle.tally_task_ = &tally_task(sr, 0);

  if (handle.is_linear_) {
    handle.centroid_ = &centroid(sr);
    handle.centroid_iteration_ = &centroid_iteration(sr);
    handle.centroid_t_ = &centroid_t(sr);
    handle.mom_matrix_ = &mom_matrix(sr);
    handle.mom_matrix_t_ = &mom_matrix_t(sr);
    handle.source_gradients_ = &source_gradients(sr, 0);
    handle.flux_moments_old_ = &flux_moments_old(sr, 0);
    handle.flux_moments_new_ = &flux_moments_new(sr, 0);
    handle.flux_moments_t_ = &flux_moments_t(sr, 0);
  }

  return handle;
}

void SourceRegionContainer::adjoint_reset()
{
  std::fill(n_hits_.begin(), n_hits_.end(), 0);
  std::fill(volume_.begin(), volume_.end(), 0.0);
  std::fill(volume_t_.begin(), volume_t_.end(), 0.0);
  std::fill(volume_sq_.begin(), volume_sq_.end(), 0.0);
  std::fill(volume_sq_t_.begin(), volume_sq_t_.end(), 0.0);
  std::fill(volume_naive_.begin(), volume_naive_.end(), 0.0);
  std::fill(
    external_source_present_.begin(), external_source_present_.end(), 0);
  std::fill(external_source_.begin(), external_source_.end(), 0.0);
  std::fill(centroid_.begin(), centroid_.end(), Position {0.0, 0.0, 0.0});
  std::fill(centroid_iteration_.begin(), centroid_iteration_.end(),
    Position {0.0, 0.0, 0.0});
  std::fill(centroid_t_.begin(), centroid_t_.end(), Position {0.0, 0.0, 0.0});
  std::fill(mom_matrix_.begin(), mom_matrix_.end(),
    MomentMatrix {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  std::fill(mom_matrix_t_.begin(), mom_matrix_t_.end(),
    MomentMatrix {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  for (auto& task_set : volume_task_) {
    task_set.clear();
  }
  std::fill(scalar_flux_old_.begin(), scalar_flux_old_.end(), 0.0);
  std::fill(scalar_flux_new_.begin(), scalar_flux_new_.end(), 0.0);
  std::fill(scalar_flux_final_.begin(), scalar_flux_final_.end(), 0.0);
  std::fill(source_.begin(), source_.end(), 0.0f);
  std::fill(external_source_.begin(), external_source_.end(), 0.0f);

  std::fill(source_gradients_.begin(), source_gradients_.end(),
    MomentArray {0.0, 0.0, 0.0});
  std::fill(flux_moments_old_.begin(), flux_moments_old_.end(),
    MomentArray {0.0, 0.0, 0.0});
  std::fill(flux_moments_new_.begin(), flux_moments_new_.end(),
    MomentArray {0.0, 0.0, 0.0});
  std::fill(flux_moments_t_.begin(), flux_moments_t_.end(),
    MomentArray {0.0, 0.0, 0.0});

  for (auto& task_set : tally_task_) {
    task_set.clear();
  }
}

} // namespace openmc
