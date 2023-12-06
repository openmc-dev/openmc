#include "openmc/random_ray/iteration.h"
#include "openmc/random_ray/tally_convert.h"
#include "openmc/random_ray/vtk_plot.h"
#include "openmc/output.h"
#include "openmc/geometry.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/simulation.h"
#include "openmc/source.h"
#include "openmc/eigenvalue.h"
#include "openmc/timer.h"
#include "openmc/mgxs_interface.h"
#include "openmc/message_passing.h"
#include "openmc/plot.h"

namespace openmc {

void all_reduce_random_ray_batch_results(bool mapped_all_tallies)
{
#ifdef OPENMC_MPI

  // If we only have 1 MPI rank, no need
  // to reduce anything.
  if (mpi::n_procs <= 1)
    return;
    
  simulation::time_bank_sendrecv.start();
  
  // The "position_recorded" variable needs to be allreduced (and maxed),
  // as whether or not a cell was hit will affect some decisions in how the
  // source is calculated in the next iteration so as to avoid dividing
  // by zero. We take the max rather than the sum as the hit values are
  // expected to be zero or 1.
  MPI_Allreduce(
      MPI_IN_PLACE,
      random_ray::position_recorded.data(),
      random_ray::n_source_regions,
      MPI_INT,
      MPI_MAX,
      mpi::intracomm);

  // The posiiton variable is more complicated to reduce than the others,
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
  
  // First, we broadcast the fully mapped tally status variable so that
  // all ranks are on the same page
  int mapped_all_tallies_i = static_cast<int>(mapped_all_tallies);
  MPI_Bcast(
      &mapped_all_tallies_i,
      1,
      MPI_INT,
      0, 
      mpi::intracomm);
  
  // Then, we perform the gather of position data, if needed
  if (simulation::current_batch > settings::n_inactive && !mapped_all_tallies_i) {
    
    // Master rank will gather results and pick valid positions
    if (mpi::master) {
      // Initialize temporary vector for receiving positions
      std::vector<std::vector<Position>> all_position;
      all_position.resize(mpi::n_procs);
      for (int i = 0; i < mpi::n_procs; i++) {
        all_position[i].resize(random_ray::n_source_regions);
      }

      // Copy master rank data into gathered vector for convenience
      all_position[0] = random_ray::position;

      // Receive all data into gather vector
      for (int i = 1; i < mpi::n_procs; i++) {
        MPI_Recv(all_position[i].data(), 
            random_ray::n_source_regions * 3,
            MPI_DOUBLE,
            i,
            0,
            mpi::intracomm,
            MPI_STATUS_IGNORE);
      }
      
      // Scan through gathered data and pick first valid cell posiiton
      for (int sr = 0; sr < random_ray::n_source_regions; sr++) {
        if (random_ray::position_recorded[sr] == 1) {
          for (int i = 0; i < mpi::n_procs; i++) {
            if (all_position[i][sr].x != 0.0 || all_position[i][sr].y != 0.0 || all_position[i][sr].z != 0.0) {
              random_ray::position[sr] = all_position[i][sr];
              break;
            }
          }
        }
      }
    } else {
      // Other ranks just send in their data
      MPI_Send(random_ray::position.data(),
          random_ray::n_source_regions * 3,
          MPI_DOUBLE,
          0,
          0,
          mpi::intracomm);
    }
  }

  // For the rest of the source region data, we simply perform an all reduce,
  // as these values will be needed on all ranks for transport during the
  // next iteration.
  MPI_Allreduce(
      MPI_IN_PLACE,
      random_ray::volume.data(),
      random_ray::n_source_regions,
      MPI_DOUBLE,
      MPI_SUM,
      mpi::intracomm);

  MPI_Allreduce(
      MPI_IN_PLACE,
      random_ray::was_hit.data(),
      random_ray::n_source_regions,
      MPI_INT,
      MPI_SUM,
      mpi::intracomm);

  MPI_Allreduce(
      MPI_IN_PLACE,
      random_ray::scalar_flux_new.data(),
      random_ray::n_source_elements,
      MPI_FLOAT,
      MPI_SUM,
      mpi::intracomm);
  
  simulation::time_bank_sendrecv.stop();
#endif
}

int openmc_run_random_ray()
{
  // Initialize OpenMC general data structures
  openmc_simulation_init();

  // Validate that inputs meet requirements for random ray mode
  if (mpi::master)
    validate_random_ray_inputs();

  // Intialize Cell (Flat Source Region) data structures
  initialize_source_regions();

  // Random ray mode does not have an inner loop over generations within a batch, so set the current gen to 1
  simulation::current_gen = 1;

  // There are no source sites in random ray mode, so be sure to disable to ensure we don't attempt
  // to write source sites to statepoint
  settings::source_write = false;

  // If all source regions have been hit by a ray during transport, then we can skip the mapping process
  // for future iterations.
  bool mapped_all_tallies = false;
  
  // Random ray eigenvalue
  // TODO: Use OpenMC native eigenvalue?
  double k_eff = 1.0;

  // Tracks the average FSR miss rate for analysis and reporting
  double avg_miss_rate = 0.0;

  // Tracks the total number of geometric intersections by all rays for reporting
  uint64_t total_geometric_intersections = 0;

  // Begin main simulation timer
  simulation::time_total.start();

  // Random ray power iteration loop
  while (simulation::current_batch < settings::n_batches) {
    
    // Initialize the current batch
    initialize_batch();
    initialize_generation();

    // Reset total starting particle weight used for normalizing tallies
    simulation::total_weight = 1.0;

    // Update source term (scattering + fission)
    update_neutron_source(k_eff);

    // Reset scalar fluxes, iteration volume tallies, and region hit flags to zero
    parallel_fill<float>(random_ray::scalar_flux_new, 0.0f);
    parallel_fill<double>(random_ray::volume, 0.0);
    parallel_fill<int>(random_ray::was_hit, 0);

    // Start timer for transport
    simulation::time_transport.start();

    // Transport sweep over all random rays for the iteration
    #pragma omp parallel for schedule(runtime) reduction(+:total_geometric_intersections)
    for (int i = 0; i < simulation::work_per_rank; i++)
    {
      RandomRay r;
      r.initialize_ray(i);
      total_geometric_intersections += r.transport_history_based_single_ray();
    }
    
    simulation::time_transport.stop();
    
    // If using multiple MPI ranks, perform all reduce on all transport results
    all_reduce_random_ray_batch_results(mapped_all_tallies);

    // Normalize scalar flux and update volumes
    normalize_scalar_flux_and_volumes();

    // Add source to scalar flux, compute number of FSR hits
    int64_t n_hits = add_source_to_scalar_flux();

    // Compute random ray k-eff
    k_eff = compute_k_eff(k_eff);

    // Store random ray k-eff into OpenMC's native k-eff variable
    global_tally_tracklength = k_eff;

    // Execute all tallying tasks, if this is an active batch
    if (simulation::current_batch > settings::n_inactive && mpi::master) {

      // Generate mapping between source regions and tallies
      if (!mapped_all_tallies) {
        mapped_all_tallies = convert_source_regions_to_tallies();
      }
      
      // Use above mapping to contribute FSR flux data to appropriate tallies
      random_ray_tally();
    
      // Add this iteration's scalar flux estimate to final accumulated estimate
      #pragma omp parallel for
      for (int64_t se = 0; se < random_ray::n_source_elements; se++) {
        random_ray::scalar_flux_final[se] += random_ray::scalar_flux_new[se];
      } 
    }
    

    // Set phi_old = phi_new
    random_ray::scalar_flux_old.swap(random_ray::scalar_flux_new);

    // Check for any obvious insabilities/nans/infs
    instability_check(n_hits, k_eff, avg_miss_rate);

    // Finalize the current batch
    finalize_generation();
    finalize_batch();
  } // End random ray power iteration loop

  openmc::simulation::time_total.stop();
  
  // Finalize OpenMC
  openmc_simulation_finalize();
  
  // Reduce number of intersections
#ifdef OPENMC_MPI
  if (mpi::n_procs > 1) {
    uint64_t total_geometric_intersections_reduced = 0;
    MPI_Reduce(
        &total_geometric_intersections,
        &total_geometric_intersections_reduced,
        1,
        MPI_UNSIGNED_LONG,
        MPI_SUM,
        0,
        mpi::intracomm);
    total_geometric_intersections = total_geometric_intersections_reduced;
  }
  #endif
  
  // Print random ray results
  if (mpi::master) {
    print_results_random_ray(total_geometric_intersections, avg_miss_rate/settings::n_batches);
    if (model::plots.size() > 0) {
      plot_3D_vtk();
    }
  }

  return 0;
}

// Compute new estimate of scattering + fission sources in each source region
// based on the flux estimate from the previous iteration.
void update_neutron_source(double k_eff)
{
  simulation::time_update_src.start();

  double inverse_k_eff = 1.0 / k_eff;
  int negroups = data::mg.num_energy_groups_;

  // Temperature and angle indices, if using multiple temperature
  // data sets and/or anisotropic data sets.
  // TODO: Currently assumes we are only using single temp/single
  // angle data.
  const int t = 0;
  const int a = 0;

  #pragma omp parallel for
  for (int sr = 0; sr < random_ray::n_source_regions; sr++) {
    int material = random_ray::material[sr]; 

    for (int energy_group_out = 0; energy_group_out < negroups; energy_group_out++) {
      float Sigma_t = data::mg.macro_xs_[material].get_xs(MgxsType::TOTAL, energy_group_out, nullptr, nullptr, nullptr, t, a);
      float scatter_source = 0.0f;
      float fission_source = 0.0f;

      for (int energy_group_in = 0; energy_group_in < negroups; energy_group_in++) {
        float scalar_flux = random_ray::scalar_flux_old[sr * negroups + energy_group_in];
        float Sigma_s = data::mg.macro_xs_[material].get_xs(MgxsType::NU_SCATTER, energy_group_in, &energy_group_out, nullptr, nullptr, t, a);
        float nu_Sigma_f = data::mg.macro_xs_[material].get_xs(MgxsType::NU_FISSION, energy_group_in, nullptr, nullptr, nullptr, t, a);
        float Chi = data::mg.macro_xs_[material].get_xs(MgxsType::CHI_PROMPT, energy_group_in, &energy_group_out, nullptr, nullptr, t, a);
        scatter_source += Sigma_s    * scalar_flux;
        fission_source += nu_Sigma_f * scalar_flux * Chi;
      }

      fission_source *= inverse_k_eff;
      float new_isotropic_source = (scatter_source + fission_source)  / Sigma_t;
      random_ray::source[sr * negroups + energy_group_out] = new_isotropic_source;
    }

  }

  simulation::time_update_src.stop();
}

// Normalizes flux and updates simulation-averaged volume estimate
void normalize_scalar_flux_and_volumes()
{
  int negroups = data::mg.num_energy_groups_;
  double total_active_distance_per_iteration = settings::random_ray_distance_active * settings::n_particles;

  float  normalization_factor =        1.0 /  total_active_distance_per_iteration;
  double volume_normalization_factor = 1.0 / (total_active_distance_per_iteration * simulation::current_batch);

  // Normalize Scalar flux to total distance travelled by all rays this iteration
  #pragma omp parallel for
  for (int64_t e = 0; e < random_ray::scalar_flux_new.size(); e++) {
    random_ray::scalar_flux_new[e] *= normalization_factor;
  }

  // Accumulate cell-wise ray length tallies collected this iteration, then
  // update the simulation-averaged cell-wise volume estimates
  #pragma omp parallel for
  for (int64_t sr = 0; sr < random_ray::n_source_regions; sr++) {
    random_ray::volume_t[sr] += random_ray::volume[sr];
    random_ray::volume[sr] = random_ray::volume_t[sr] * volume_normalization_factor;
  }
}

// Combines transport flux contributions and flat source contributions
// from the previous iteration to generate this iteration's estimate of
// scalar flux.
int64_t add_source_to_scalar_flux()
{
  int negroups = data::mg.num_energy_groups_;

  int64_t n_hits = 0;
  
  // Temperature and angle indices, if using multiple temperature
  // data sets and/or anisotropic data sets.
  // TODO: Currently assumes we are only using single temp/single
  // angle data.
  const int t = 0;
  const int a = 0;

  #pragma omp parallel for reduction(+:n_hits)
  for (int sr = 0; sr < random_ray::n_source_regions; sr++) {

    // Check if this cell was hit this iteration
    int was_cell_hit = random_ray::was_hit[sr];
    if (was_cell_hit) {
      n_hits++;
    }

    double volume = random_ray::volume[sr];
    int material = random_ray::material[sr]; 
    for (int e = 0; e < negroups; e++) {
      int64_t idx = (sr * negroups) + e;

      // There are three scenarios we need to consider:
      if (was_cell_hit) {
        // 1. If the FSR was hit this iteration, then the new flux is equal to the flat source from the previous iteration
        // plus the contributions from rays passing through the source region (computed during the transport sweep)
        float Sigma_t = data::mg.macro_xs_[material].get_xs(MgxsType::TOTAL, e, nullptr, nullptr, nullptr, t, a);
        random_ray::scalar_flux_new[idx] /= (Sigma_t * volume);
        random_ray::scalar_flux_new[idx] += random_ray::source[idx];
      } else if (volume > 0.0) {
        // 2. If the FSR was not hit this iteration, but has been hit some previous iteration, then we simply set the
        // new scalar flux to be equal to the contribution from the flat source alone. 
        random_ray::scalar_flux_new[idx] = random_ray::source[idx];
      } else {
        // If the FSR was not hit this iteration, and it has never been hit in any iteration (i.e., volume is zero),
        // then we want to set this to 0 to avoid dividing anything by a zero volume.
        random_ray::scalar_flux_new[idx] = 0.f;
      }
    }
  }

  // Return the number of source regions that were hit this iteration
  return n_hits;
}

// Generates new estimate of k_eff based on the differences between this
// iteration's estimate of the scalar flux and the last iteration's estimate.
double compute_k_eff(double k_eff_old)
{
  int negroups = data::mg.num_energy_groups_;
  double fission_rate_old = 0;
  double fission_rate_new = 0;
  
  // Temperature and angle indices, if using multiple temperature
  // data sets and/or anisotropic data sets.
  // TODO: Currently assumes we are only using single temp/single
  // angle data.
  const int t = 0;
  const int a = 0;

  #pragma omp parallel for reduction(+:fission_rate_old, fission_rate_new)
  for (int sr = 0; sr < random_ray::n_source_regions; sr++) {

    // If simulation averaged volume is zero, don't include this cell
    double volume = random_ray::volume[sr];
    if (volume == 0.0) {
      continue;
    }

    int material = random_ray::material[sr]; 

    double sr_fission_source_old = 0;
    double sr_fission_source_new = 0;

    for (int e = 0; e < negroups; e++) {
      int64_t idx = (sr * negroups) + e;
      double nu_Sigma_f = data::mg.macro_xs_[material].get_xs(MgxsType::NU_FISSION, e, nullptr, nullptr, nullptr, t, a);
      sr_fission_source_old += nu_Sigma_f * random_ray::scalar_flux_old[idx];
      sr_fission_source_new += nu_Sigma_f * random_ray::scalar_flux_new[idx];
    }

    fission_rate_old += sr_fission_source_old * volume;
    fission_rate_new += sr_fission_source_new * volume;
  }


  double k_eff_new = k_eff_old * (fission_rate_new / fission_rate_old);

  return k_eff_new;
}

// Apply a few sanity checks to catch obvious cases of numerical instability.
// Instability typically only occurs if ray density is extremely low.
void instability_check(int64_t n_hits, double k_eff, double& avg_miss_rate)
{
  double percent_missed = ((random_ray::n_source_regions - n_hits) / static_cast<double>(random_ray::n_source_regions)) * 100.0;
  avg_miss_rate += percent_missed;

  if( percent_missed > 10.0 ) {
    warning(fmt::format("Very high FSR miss rate detected ({:.3f}%). Instability may occur. Increase ray density by adding more rays and/or active distance.", percent_missed));
  } else if( percent_missed > 0.01 ) {
    warning(fmt::format("Elevated FSR miss rate detected ({:.3f}%). Increasing ray density by adding more rays and/or active distance may improve simulation efficiency.", percent_missed));
  }

  if (k_eff > 10.0 || k_eff < 0.01 || !(std::isfinite(k_eff))) {
    fatal_error("Instability detected");
  }
}

// Enforces restrictions on inputs in random ray mode.  While there are
// many features that don't make sense in random ray mode, and are therefore
// unsupported, we limit our testing/enforcement operations only to inputs
// that may cause erroneous/misleading output or crashes from the solver.
void validate_random_ray_inputs()
{
  // Validate tallies
  for (auto& tally : model::tallies) {

    // Validate score types
    for (auto score_bin : tally->scores_) {
      switch (score_bin) {
        case SCORE_FLUX:
        case SCORE_TOTAL:
        case SCORE_FISSION:
        case SCORE_NU_FISSION:
        case SCORE_EVENTS:
          break;
        default:
          fatal_error("Invalid score specified. Only flux, total, fission, nu-fission, and event scores are supported in random ray mode.");
      }
    }

    // Validate estimator types
    if (tally->estimator_ != TallyEstimator::ANALOG) {
      fatal_error("Invalid estimator specified. Only analog estimators are supported in random ray mode.");
    }

    // Validate filter types
    for (auto f : tally->filters()) {
      auto& filter = *model::tally_filters[f];

      switch (filter.type()) {
        case  FilterType::CELL:
        case  FilterType::CELL_INSTANCE:
        case  FilterType::DISTRIBCELL:
        case  FilterType::ENERGY:
        case  FilterType::MATERIAL:
        case  FilterType::MESH:
        case  FilterType::UNIVERSE:
          break;
        default:
          fatal_error("Invalid filter specified. Only cell, cell_instance, distribcell, energy, material, mesh, and universe filters are supported in random ray mode.");
      }
    }
  }

  // Validate MGXS data
  for (auto& material : data::mg.macro_xs_) {
    if (!material.is_isotropic) {
      fatal_error("Anisotropic MGXS detected. Only isotropic XS data sets supported in random ray mode.");
    }
    if (material.get_xsdata().size() > 1) {
      fatal_error("Non-isothermal MGXS detected. Only isothermal XS data sets supported in random ray mode.");
    }
  }

  // Validate solver mode
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    fatal_error("Invalid run mode. Fixed source not yet supported in random ray mode.");
  }
  
  // Validate eigenvalue sources (must be single uniform isotropic box source)

  // Ensure only a single source is present
  if (model::external_sources.size() != 1) {
    fatal_error("Invalid source definition -- only a single source is allowed in random ray mode.");
  }

  Source* s = model::external_sources[0].get();

  // Check for independent source
  IndependentSource* is = dynamic_cast<IndependentSource*>(s);
  if (is == nullptr) {
    fatal_error("Invalid source definition -- only independent sources are allowed in random ray mode.");
  }

  // Check for box source
  SpatialDistribution* space_dist = is->space();
  SpatialBox* sb = dynamic_cast<SpatialBox*>(space_dist);
  if (sb == nullptr) {
    fatal_error("Invalid source definition -- only box sources are allowed in random ray mode. If no source is specified, OpenMC default is an isotropic point source at the origin, which is invalid in random ray mode.");
  }

  // Check that box source is not restricted to fissionable areas
  if (sb->only_fissionable()) {
    fatal_error("Invalid source definition -- fissionable spatial distribution not allowed in random ray mode.");
  }

  // Check for isotropic source
  UnitSphereDistribution* angle_dist = is->angle();
  Isotropic* id = dynamic_cast<Isotropic*>(angle_dist);
  if (id == nullptr) {
    fatal_error("Invalid source definition -- only isotropic sources are allowed in random ray mode.");
  }

  // Check plotting files
  for (int p = 0; p < model::plots.size(); p++) {

    // Get handle to OpenMC plot object and extract params
    Plot* openmc_plot = dynamic_cast<Plot*>(model::plots[p].get());
      
    // Random ray plots only support voxel plots
    if (openmc_plot == nullptr) {
      warning(fmt::format("Plot {} will not be used for end of simulation data plotting -- only voxel plotting is allowed in random ray mode.", p));
      continue;
    } else if (openmc_plot->type_ != Plot::PlotType::voxel) {
      warning(fmt::format("Plot {} will not be used for end of simulation data plotting -- only voxel plotting is allowed in random ray mode.", p));
      continue;
    } 
  }

  // Warn about slow MPI domain replication, if detected
  #ifdef OPENMC_MPI
  if (mpi::n_procs > 1) {
    warning("Domain replication in random ray is supported, but suffers from poor scaling of source all-reduce operations. Performance may severely degrade beyond just a few MPI ranks. Domain decomposition may be implemented in the future to provide efficient scaling.");
  }
  #endif
}

} // namespace openmc
