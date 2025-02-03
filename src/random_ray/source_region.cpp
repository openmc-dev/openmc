#include "openmc/random_ray/source_region.h"

#include "openmc/error.h"
#include "openmc/message_passing.h"

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

SourceRegion::SourceRegion(const SourceRegionHandle& handle)
  : SourceRegion(handle.negroups_, handle.mom_matrix_ != nullptr)
{
  material_ = handle.material();
  external_source_present_ = handle.external_source_present();
  mesh_ = handle.mesh();
  for (int g = 0; g < scalar_flux_new_.size(); g++) {
    scalar_flux_old_[g] = handle.scalar_flux_old_[g];
    source_[g] = handle.source_[g];
    external_source_[g] = handle.external_source_[g];
  }
}

SourceRegionHandle SourceRegion::get_source_region_handle()
{
  SourceRegionHandle handle;
  handle.negroups_ = scalar_flux_old_.size();
  handle.material_ = &material_;
  handle.lock_ = &lock_;
  handle.volume_ = &volume_;
  handle.volume_t_ = &volume_t_;
  handle.volume_naive_ = &volume_naive_;
  handle.position_recorded_ = &position_recorded_;
  handle.external_source_present_ = &external_source_present_;
  handle.position_ = &position_;
  handle.centroid_ = &centroid_;
  handle.centroid_iteration_ = &centroid_iteration_;
  handle.centroid_t_ = &centroid_t_;
  handle.mom_matrix_ = &mom_matrix_;
  handle.mom_matrix_t_ = &mom_matrix_t_;
  handle.volume_task_ = &volume_task_;
  handle.mesh_ = &mesh_;
  handle.scalar_flux_old_ = scalar_flux_old_.data();
  handle.scalar_flux_new_ = scalar_flux_new_.data();
  handle.source_ = source_.data();
  handle.external_source_ = external_source_.data();
  handle.scalar_flux_final_ = scalar_flux_final_.data();
  handle.source_gradients_ = source_gradients_.data();
  handle.flux_moments_old_ = flux_moments_old_.data();
  handle.flux_moments_new_ = flux_moments_new_.data();
  handle.flux_moments_t_ = flux_moments_t_.data();
  handle.tally_task_ = tally_task_.data();
  return handle;
}

//==============================================================================
// SourceRegionContainer implementation
//==============================================================================

void SourceRegionContainer::push_back(const SourceRegion& sr)
{
  n_source_regions_++;

  // Scalar fields
  material_.push_back(sr.material_);
  lock_.push_back(sr.lock_);
  volume_.push_back(sr.volume_);
  volume_t_.push_back(sr.volume_t_);
  volume_naive_.push_back(sr.volume_naive_);
  position_recorded_.push_back(sr.position_recorded_);
  external_source_present_.push_back(sr.external_source_present_);
  position_.push_back(sr.position_);
  volume_task_.push_back(sr.volume_task_);
  mesh_.push_back(sr.mesh_);

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
  lock_.clear();
  volume_.clear();
  volume_t_.clear();
  volume_naive_.clear();
  position_recorded_.clear();
  external_source_present_.clear();
  position_.clear();
  mesh_.clear();

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

void SourceRegionContainer::mpi_sync_ranks(bool reduce_position)
{
#ifdef OPENMC_MPI

  // The "position_recorded" variable needs to be allreduced (and maxed),
  // as whether or not a cell was hit will affect some decisions in how the
  // source is calculated in the next iteration so as to avoid dividing
  // by zero. We take the max rather than the sum as the hit values are
  // expected to be zero or 1.
  MPI_Allreduce(MPI_IN_PLACE, position_recorded_.data(), 1, MPI_INT, MPI_MAX,
    mpi::intracomm);

  // The position variable is more complicated to reduce than the others,
  // as we do not want the sum of all positions in each cell, rather, we
  // want to just pick any single valid position. Thus, we perform a gather
  // and then pick the first valid position we find for all source regions
  // that have had a position recorded. This operation does not need to
  // be broadcast back to other ranks, as this value is only used for the
  // tally conversion operation, which is only performed on the master rank.
  // While this is expensive, it only needs to be done for active batches,
  // and only if we have not mapped all the tallies yet. Once tallies are
  // fully mapped, then the position vector is fully populated, so this
  // operation can be skipped.

  // Then, we perform the gather of position data, if needed
  if (reduce_position) {

    // Master rank will gather results and pick valid positions
    if (mpi::master) {
      // Initialize temporary vector for receiving positions
      vector<vector<Position>> all_position;
      all_position.resize(mpi::n_procs);
      for (int i = 0; i < mpi::n_procs; i++) {
        all_position[i].resize(n_source_regions_);
      }

      // Copy master rank data into gathered vector for convenience
      all_position[0] = position_;

      // Receive all data into gather vector
      for (int i = 1; i < mpi::n_procs; i++) {
        MPI_Recv(all_position[i].data(), n_source_regions_ * 3, MPI_DOUBLE, i,
          0, mpi::intracomm, MPI_STATUS_IGNORE);
      }

      // Scan through gathered data and pick first valid cell posiiton
      for (int64_t sr = 0; sr < n_source_regions_; sr++) {
        if (position_recorded_[sr] == 1) {
          for (int i = 0; i < mpi::n_procs; i++) {
            if (all_position[i][sr].x != 0.0 || all_position[i][sr].y != 0.0 ||
                all_position[i][sr].z != 0.0) {
              position_[sr] = all_position[i][sr];
              break;
            }
          }
        }
      }
    } else {
      // Other ranks just send in their data
      MPI_Send(position_.data(), n_source_regions_ * 3, MPI_DOUBLE, 0, 0,
        mpi::intracomm);
    }
  }

  // For the rest of the source region data, we simply perform an all reduce,
  // as these values will be needed on all ranks for transport during the
  // next iteration.
  MPI_Allreduce(MPI_IN_PLACE, volume_.data(), n_source_regions_, MPI_DOUBLE,
    MPI_SUM, mpi::intracomm);

  MPI_Allreduce(MPI_IN_PLACE, scalar_flux_new_.data(),
    n_source_regions_ * negroups_, MPI_DOUBLE, MPI_SUM, mpi::intracomm);

  if (is_linear_) {
    // We are going to assume we can safely cast Position, MomentArray,
    // and MomentMatrix to contiguous arrays of doubles for the MPI
    // allreduce operation. This is a safe assumption as typically
    // compilers will at most pad to 8 byte boundaries. If a new FP32
    // MomentArray type is introduced, then there will likely be padding, in
    // which case this function will need to become more complex.
    if (sizeof(MomentArray) != 3 * sizeof(double) ||
        sizeof(MomentMatrix) != 6 * sizeof(double)) {
      fatal_error(
        "Unexpected buffer padding in linear source domain reduction.");
    }

    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(flux_moments_new_.data()),
      n_source_regions_ * negroups_ * 3, MPI_DOUBLE, MPI_SUM, mpi::intracomm);
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(mom_matrix_.data()),
      n_source_regions_ * 6, MPI_DOUBLE, MPI_SUM, mpi::intracomm);
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(centroid_iteration_.data()),
      n_source_regions_ * 3, MPI_DOUBLE, MPI_SUM, mpi::intracomm);
  }

#endif
}

void SourceRegionContainer::reduce_to_base()
{
  lock_.clear();
  volume_.clear();
  volume_t_.clear();
  volume_naive_.clear();
  position_recorded_.clear();
  position_.clear();
  centroid_.clear();
  centroid_iteration_.clear();
  centroid_t_.clear();
  mom_matrix_.clear();
  mom_matrix_t_.clear();
  volume_task_.clear();
  scalar_flux_new_.clear();
  scalar_flux_final_.clear();
  source_gradients_.clear();
  flux_moments_old_.clear();
  flux_moments_new_.clear();
  flux_moments_t_.clear();
  tally_task_.clear();
}

SourceRegionHandle SourceRegionContainer::get_source_region_handle(int64_t sr)
{
  SourceRegionHandle handle;
  handle.negroups_ = negroups();
  handle.material_ = &material(sr);
  handle.lock_ = &lock(sr);
  handle.volume_ = &volume(sr);
  handle.volume_t_ = &volume_t(sr);
  handle.volume_naive_ = &volume_naive(sr);
  handle.position_recorded_ = &position_recorded(sr);
  handle.external_source_present_ = &external_source_present(sr);
  handle.position_ = &position(sr);
  handle.volume_task_ = &volume_task(sr);
  handle.mesh_ = &mesh(sr);
  handle.scalar_flux_old_ = &scalar_flux_old(sr, 0);
  handle.scalar_flux_new_ = &scalar_flux_new(sr, 0);
  handle.source_ = &source(sr, 0);
  handle.external_source_ = &external_source(sr, 0);
  handle.scalar_flux_final_ = &scalar_flux_final(sr, 0);
  handle.tally_task_ = &tally_task(sr, 0);

  if (is_linear_) {
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

} // namespace openmc