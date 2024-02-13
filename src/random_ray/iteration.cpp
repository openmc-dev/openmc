#include "openmc/random_ray/iteration.h"
#include "openmc/eigenvalue.h"
#include "openmc/geometry.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/output.h"
#include "openmc/plot.h"
#include "openmc/random_ray/tally_convert.h"
#include "openmc/random_ray/vtk_plot.h"
#include "openmc/simulation.h"
#include "openmc/source.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/timer.h"

namespace openmc {


RandomRaySimulation::RandomRaySimulation()
{
  // Determine which source term should be used for ray sampling
  find_random_ray_sampling_source_index();
  
  // There are no source sites in random ray mode, so be sure to disable to
  // ensure we don't attempt to write source sites to statepoint
  settings::source_write = false;
}

void RandomRaySimulation::run_simulation()
{
  // Random ray power iteration loop
  while (simulation::current_batch < settings::n_batches) {

    // Initialize the current batch
    initialize_batch();
    initialize_generation();

    // Reset total starting particle weight used for normalizing tallies
    simulation::total_weight = 1.0;

    // Update source term (scattering + fission)
    domain_.update_neutron_source();

    // Reset scalar fluxes, iteration volume tallies, and region hit flags to
    // zero
    domain_.batch_reset();

    // Start timer for transport
    simulation::time_transport.start();

// Transport sweep over all random rays for the iteration
#pragma omp parallel for schedule(dynamic)                                     \
  reduction(+ : total_geometric_intersections)
    for (int i = 0; i < simulation::work_per_rank; i++) {
      RandomRay ray;
      ray.initialize_ray(i, sampling_source_, domain_);
      total_geometric_intersections += ray.transport_history_based_single_ray();
    }

    simulation::time_transport.stop();

    // If using multiple MPI ranks, perform all reduce on all transport results
    domain_.all_reduce_random_ray_batch_results();

    // Normalize scalar flux and update volumes
    domain_.normalize_scalar_flux_and_volumes();

    // Add source to scalar flux, compute number of FSR hits
    int64_t n_hits = domain_.add_source_to_scalar_flux();

    // Compute random ray k-eff
    k_eff_ = domain_.compute_k_eff(k_eff_);

    // Store random ray k-eff into OpenMC's native k-eff variable
    global_tally_tracklength = k_eff_;

    // Execute all tallying tasks, if this is an active batch
    if (simulation::current_batch > settings::n_inactive && mpi::master) {

      // Generate mapping between source regions and tallies
      if (!domain_.mapped_all_tallies) {
        domain_.convert_source_regions_to_tallies();
      }

      // Use above mapping to contribute FSR flux data to appropriate tallies
      domain_.random_ray_tally();

      // Add this iteration's scalar flux estimate to final accumulated estimate
      domain_.accumulate_iteration_flux();
    }

    // Set phi_old = phi_new
    domain_.scalar_flux_old_.swap(domain_.scalar_flux_new_);

    // Check for any obvious insabilities/nans/infs
    instability_check(n_hits, k_eff_, avg_miss_rate_);

    // Finalize the current batch
    finalize_generation();
    finalize_batch();
  } // End random ray power iteration loop
}

void RandomRaySimulation::reduce_simulation_statistics()
{
  // Reduce number of intersections
#ifdef OPENMC_MPI
  if (mpi::n_procs > 1) {
    uint64_t total_geometric_intersections_reduced = 0;
    MPI_Reduce(&total_geometric_intersections_,
      &total_geometric_intersections_reduced, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0,
      mpi::intracomm);
    total_geometric_intersections_ = total_geometric_intersections_reduced;
  }
#endif
}

void RandomRaySimulation::output_simulation_data()
{
  // Print random ray results
  if (mpi::master) {
    print_results_random_ray(
      total_geometric_intersections, avg_miss_rate / settings::n_batches);
    if (model::plots.size() > 0) {
      plot_3D_vtk();
    }
  }
}

void openmc_run_random_ray()
{
  // Initialize OpenMC general data structures
  openmc_simulation_init();

  // Validate that inputs meet requirements for random ray mode
  if (mpi::master)
    validate_random_ray_inputs();
  
  // Initialize Random Ray Simulation Object
  RandomRaySimulation sim;

  // Random ray mode does not have an inner loop over generations within a
  // batch, so set the current gen to 1
  simulation::current_gen = 1;

  // Begin main simulation timer
  simulation::time_total.start();

  // Execute random ray simulation
  sim.run_simulation();

  // End main simulation timer
  openmc::simulation::time_total.stop();

  // Finalize OpenMC
  openmc_simulation_finalize();
}


// Apply a few sanity checks to catch obvious cases of numerical instability.
// Instability typically only occurs if ray density is extremely low.
void RandomRaySimulation::instability_check(int64_t n_hits, double k_eff, double& avg_miss_rate)
{
  double percent_missed = ((domain_.n_source_regions - n_hits) /
                            static_cast<double>(domain_.n_source_regions)) *
                          100.0;
  avg_miss_rate += percent_missed;

  if (percent_missed > 10.0) {
    warning(fmt::format(
      "Very high FSR miss rate detected ({:.3f}%). Instability may occur. "
      "Increase ray density by adding more rays and/or active distance.",
      percent_missed));
  } else if (percent_missed > 0.01) {
    warning(fmt::format("Elevated FSR miss rate detected ({:.3f}%). Increasing "
                        "ray density by adding more rays and/or active "
                        "distance may improve simulation efficiency.",
      percent_missed));
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
  ///////////////////////////////////////////////////////////////////
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
        fatal_error(
          "Invalid score specified. Only flux, total, fission, nu-fission, and "
          "event scores are supported in random ray mode.");
      }
    }

    // Validate estimator types
    if (tally->estimator_ != TallyEstimator::ANALOG) {
      fatal_error("Invalid estimator specified. Only analog estimators are "
                  "supported in random ray mode.");
    }

    // Validate filter types
    for (auto f : tally->filters()) {
      auto& filter = *model::tally_filters[f];

      switch (filter.type()) {
      case FilterType::CELL:
      case FilterType::CELL_INSTANCE:
      case FilterType::DISTRIBCELL:
      case FilterType::ENERGY:
      case FilterType::MATERIAL:
      case FilterType::MESH:
      case FilterType::UNIVERSE:
        break;
      default:
        fatal_error("Invalid filter specified. Only cell, cell_instance, "
                    "distribcell, energy, material, mesh, and universe filters "
                    "are supported in random ray mode.");
      }
    }
  }

  // Validate MGXS data
  ///////////////////////////////////////////////////////////////////
  for (auto& material : data::mg.macro_xs_) {
    if (!material.is_isotropic) {
      fatal_error("Anisotropic MGXS detected. Only isotropic XS data sets "
                  "supported in random ray mode.");
    }
    if (material.get_xsdata().size() > 1) {
      fatal_error("Non-isothermal MGXS detected. Only isothermal XS data sets "
                  "supported in random ray mode.");
    }
  }

  // Validate solver mode
  ///////////////////////////////////////////////////////////////////
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    fatal_error(
      "Invalid run mode. Fixed source not yet supported in random ray mode.");
  }

  // Validate sources
  ///////////////////////////////////////////////////////////////////

  int n_random_ray_sources = 0;

  for (int i = 0; i < model::external_sources.size(); i++) {
    Source* s = model::external_sources[i].get();

    // Check for independent source
    IndependentSource* is = dynamic_cast<IndependentSource*>(s);

    // Skip source if it is not independent, as this implies it is not
    // the random ray source
    if (is == nullptr) {
      continue;
    }

    // Skip source if this is not a random ray source
    if (is->particle_type() != ParticleType::random_ray) {
      continue;
    }

    // Increment random ray source counter
    n_random_ray_sources++;

    // Check for box source
    SpatialDistribution* space_dist = is->space();
    SpatialBox* sb = dynamic_cast<SpatialBox*>(space_dist);
    if (sb == nullptr) {
      fatal_error(
        "Invalid source definition -- only box sources are allowed in random "
        "ray "
        "mode. If no source is specified, OpenMC default is an isotropic point "
        "source at the origin, which is invalid in random ray mode.");
    }

    // Check that box source is not restricted to fissionable areas
    if (sb->only_fissionable()) {
      fatal_error(
        "Invalid source definition -- fissionable spatial distribution "
        "not allowed for random ray source.");
    }

    // Check for isotropic source
    UnitSphereDistribution* angle_dist = is->angle();
    Isotropic* id = dynamic_cast<Isotropic*>(angle_dist);
    if (id == nullptr) {
      fatal_error("Invalid source definition -- only isotropic sources are "
                  "allowed for random ray source.");
    }
  }

  // Ensure that exactly one random ray source was defined
  if (n_random_ray_sources != 1) {
    fatal_error("Invalid source definition -- a single random ray source type "
                "must be provided in random ray mode.");
  }

  // Validate plotting files
  ///////////////////////////////////////////////////////////////////
  for (int p = 0; p < model::plots.size(); p++) {

    // Get handle to OpenMC plot object
    Plot* openmc_plot = dynamic_cast<Plot*>(model::plots[p].get());

    // Random ray plots only support voxel plots
    if (openmc_plot == nullptr) {
      warning(fmt::format(
        "Plot {} will not be used for end of simulation data plotting -- only "
        "voxel plotting is allowed in random ray mode.",
        p));
      continue;
    } else if (openmc_plot->type_ != Plot::PlotType::voxel) {
      warning(fmt::format(
        "Plot {} will not be used for end of simulation data plotting -- only "
        "voxel plotting is allowed in random ray mode.",
        p));
      continue;
    }
  }

  // Warn about slow MPI domain replication, if detected
  ///////////////////////////////////////////////////////////////////
#ifdef OPENMC_MPI
  if (mpi::n_procs > 1) {
    warning(
      "Domain replication in random ray is supported, but suffers from poor "
      "scaling of source all-reduce operations. Performance may severely "
      "degrade beyond just a few MPI ranks. Domain decomposition may be "
      "implemented in the future to provide efficient scaling.");
  }
#endif
}

// Determines which external source to use for sampling
// ray starting locations/directions. Assumes sources have
// already been validated.
void RandomRaySimulation::find_random_ray_sampling_source_index()
{
  int i = 0;
  for (; i < model::external_sources.size(); i++) {
    Source* s = model::external_sources[i].get();

    // Check for independent source
    IndependentSource* is = dynamic_cast<IndependentSource*>(s);

    // Skip source if it is not independent, as this implies it is not
    // the random ray source
    if (is == nullptr) {
      continue;
    }

    // Return index when random ray source is found
    if (is->particle_type() == ParticleType::random_ray) {
      break;
    }
  }
  sampling_source_ = i;;
}

} // namespace openmc
