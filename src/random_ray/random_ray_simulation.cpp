#include "openmc/random_ray/random_ray_simulation.h"

#include "openmc/eigenvalue.h"
#include "openmc/geometry.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/output.h"
#include "openmc/plot.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/random_ray.h"
#include "openmc/simulation.h"
#include "openmc/source.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/timer.h"

namespace openmc {

//==============================================================================
// Non-member functions
//==============================================================================

void openmc_run_random_ray()
{
  // Initialize OpenMC general data structures
  openmc_simulation_init();

  // Validate that inputs meet requirements for random ray mode
  if (mpi::master)
    validate_random_ray_inputs();

  // Initialize Random Ray Simulation Object
  RandomRaySimulation sim;

  // Begin main simulation timer
  simulation::time_total.start();

  // Execute random ray simulation
  sim.simulate();

  // End main simulation timer
  openmc::simulation::time_total.stop();

  // Finalize OpenMC
  openmc_simulation_finalize();

  // Reduce variables across MPI ranks
  sim.reduce_simulation_statistics();

  // Output all simulation results
  sim.output_simulation_results();
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

  // Validate ray source
  ///////////////////////////////////////////////////////////////////

  // Check for independent source
  IndependentSource* is =
    dynamic_cast<IndependentSource*>(RandomRay::ray_source_.get());
  if (!is) {
    fatal_error("Invalid ray source definition. Ray source must provided and "
                "be of type IndependentSource.");
  }

  // Check for box source
  SpatialDistribution* space_dist = is->space();
  SpatialBox* sb = dynamic_cast<SpatialBox*>(space_dist);
  if (!sb) {
    fatal_error(
      "Invalid ray source definition -- only box sources are allowed.");
  }

  // Check that box source is not restricted to fissionable areas
  if (sb->only_fissionable()) {
    fatal_error(
      "Invalid ray source definition -- fissionable spatial distribution "
      "not allowed.");
  }

  // Check for isotropic source
  UnitSphereDistribution* angle_dist = is->angle();
  Isotropic* id = dynamic_cast<Isotropic*>(angle_dist);
  if (!id) {
    fatal_error("Invalid ray source definition -- only isotropic sources are "
                "allowed.");
  }

  // Validate external sources
  ///////////////////////////////////////////////////////////////////
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    if (model::external_sources.size() < 1) {
      fatal_error("Must provide a particle source (in addition to ray source) "
                  "in fixed source random ray mode.");
    }

    for (int i = 0; i < model::external_sources.size(); i++) {
      Source* s = model::external_sources[i].get();

      // Check for independent source
      IndependentSource* is = dynamic_cast<IndependentSource*>(s);

      if (!is) {
        fatal_error(
          "Only IndependentSource external source types are allowed in "
          "random ray mode");
      }

      // Check for isotropic source
      UnitSphereDistribution* angle_dist = is->angle();
      Isotropic* id = dynamic_cast<Isotropic*>(angle_dist);
      if (!id) {
        fatal_error(
          "Invalid source definition -- only isotropic external sources are "
          "allowed in random ray mode.");
      }

      // Validate that a domain ID was specified
      if (is->domain_ids().size() == 0) {
        if (!settings::FIRST_COLLIDED_FLUX) {
          fatal_error("Fixed sources must be specified by domain "
                      "id (cell, material, or universe) in random ray mode.");
        } else if (settings::FIRST_COLLIDED_FLUX) {
          //
          fmt::print("??? add First collided condition and fatal error here\n");
        }
      }

      // Check that a discrete energy distribution was used
      Distribution* d = is->energy();
      Discrete* dd = dynamic_cast<Discrete*>(d);
      if (!dd) {
        fatal_error(
          "Only discrete (multigroup) energy distributions are allowed for "
          "external sources in random ray mode.");
      }
    }
  }

  // Validate plotting files
  ///////////////////////////////////////////////////////////////////
  for (int p = 0; p < model::plots.size(); p++) {

    // Get handle to OpenMC plot object
    Plot* openmc_plot = dynamic_cast<Plot*>(model::plots[p].get());

    // Random ray plots only support voxel plots
    if (!openmc_plot) {
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

//==============================================================================
// RandomRaySimulation implementation
//==============================================================================

RandomRaySimulation::RandomRaySimulation()
  : negroups_(data::mg.num_energy_groups_)
{
  // There are no source sites in random ray mode, so be sure to disable to
  // ensure we don't attempt to write source sites to statepoint
  settings::source_write = false;

  // Random ray mode does not have an inner loop over generations within a
  // batch, so set the current gen to 1
  simulation::current_gen = 1;

  switch (RandomRay::source_shape_) {
  case RandomRaySourceShape::FLAT:
    domain_ = make_unique<FlatSourceDomain>();
    break;
  case RandomRaySourceShape::LINEAR:
    domain_ = make_unique<LinearSourceDomain>();
    break;
  case RandomRaySourceShape::LINEAR_XY:
    domain_ = make_unique<LinearSourceDomain>();
    break;
  default:
    fatal_error("Unknown random ray source shape");
  }
}

void RandomRaySimulation::simulate()
{
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    if (settings::FIRST_COLLIDED_FLUX == false) {
      // Transfer external source user inputs onto random ray source regions
      domain_->convert_external_sources();
      domain_->count_external_source_regions();
    } else {
      // FIRST_COLLIDED_FLUX calculation routine:
      fmt::print("==============>     FIRST COLLIDED SOURCE CONVERSION     "
                 "<=================\n\n");

      //==============================================================================
      // Turn on volume pre calculation
      double start_volume_estimation = omp_get_wtime();
      fmt::print("Computing starting volume estimates... \n");
      RandomRay::uncollided_flux_volume = true;
      simulation::current_batch = 1;

      // Pre calculate volume
#pragma omp parallel for schedule(dynamic)                                     \
  reduction(+ : total_geometric_intersections_)
      for (int i = 0; i < settings::n_volume_estimator_rays; i++) {
        RandomRay ray(i, domain_.get());
        total_geometric_intersections_ +=
          ray.transport_history_based_single_ray();
      }

      // fmt::print("simulation_batch = {:}", simulation::current_batch);
      //  Normalizing volumes
      domain_->normalize_scalar_flux_and_volumes(
        settings::n_volume_estimator_rays * RandomRay::distance_active_);

      // Print volume estimation simulation time
      double end_volume_estimation = omp_get_wtime();
      fmt::print("Volume estimation completed.\n");
      fmt::print("Run time = {:.4f} \n\n",
        end_volume_estimation - start_volume_estimation);

      // Reset parameters for First Collided Method
      RandomRay::uncollided_flux_volume = false;
      domain_->batch_reset();

      //==============================================================================
      // First Collided Transport Loop
      // It will operate until (1-3):
      // (1) There is no new FSR hits from one (n-1) iteration to the next (n)
      // (2) Reached pre-set maximum n_uncollided_rays
      // (3) Hit 100% of the FSRs
      fmt::print("  Batch    Rays         Total Source Regions Discovered\n"
                 "  ======   ==========   ================================\n");
      while (domain_->new_fsr_fc) {

        // Ray tracing and attenuation
#pragma omp parallel for schedule(dynamic)                                     \
  reduction(+ : total_geometric_intersections_)
        for (int i = old_n_rays; i < new_n_rays; i++) {
          RandomRay ray(i, domain_.get());
          total_geometric_intersections_ +=
            ray.transport_history_based_single_ray_first_collided();
        }

        int64_t n_hits_new = domain_->check_fsr_hits();

        // print results
        fmt::print("  {:6}   {:10}   {:}\n", batch_first_collided, new_n_rays,
          n_hits_new);

        // BREAK STATEMENT if Max rays reached or 100% FSR HIT

        if (n_hits_new == n_hits_old) {
          domain_->new_fsr_fc = {false};
        } else if (new_n_rays >= n_rays_max) {
          domain_->new_fsr_fc = {false};
        } else if (static_cast<double>(n_hits_new) >=
                   domain_->n_source_regions_) {
          domain_->new_fsr_fc = {false};
        }
        // Section to add in n_source_calculations
        // domain_->normalize_uncollided_scalar_flux(new_n_rays);
        // domain_->uncollided_sum_source();

        // update values for next batch
        old_n_rays = new_n_rays;
        new_n_rays *= 2;
        batch_first_collided++;
        n_hits_old = n_hits_new;
      }

      // fmt::print("old_n_rays = {:}    and new_n_rays = {:}\n", old_n_rays,
      // new_n_rays);
      //  Compute scalar_uncollided_flux
      domain_->compute_uncollided_scalar_flux();

      // Normalize scalar_uncollided_flux
      domain_->normalize_uncollided_scalar_flux(old_n_rays);

      // compute first_collided_fixed_source -
      domain_->compute_first_collided_flux();

      // Compute external source a priori.
      if (settings::volume_online_option == 1) {
        domain_->update_external_linear_source();
        domain_->update_external_flat_source();   // comment if option 2
        domain_->update_volume_uncollided_flux(); // comment if option 2
        fmt::print("Flux and moments are calculated a priori to RR \n");
      } else if (settings::volume_online_option == 2) {
        domain_->update_external_linear_source();
        fmt::print("Moments are calculated a priori to RR \n");
      } else {
        fmt::print("Flux and moments are updated online in RR \n");
      }

      // Reset Terms to RR method
      domain_->batch_reset(); // clean-up of key variables (preserves volume)
      simulation::current_batch = 0; // garantee the first batch will be 1 in RR

      // reset values for RandomRay iteration
      domain_->first_collided_mode = {true};   // add FC contribution to RR
      settings::FIRST_COLLIDED_FLUX = {false}; // regular RR calculations

      // fmt::print("first_collided_mode = {:}\n",
      // domain_->first_collided_mode);
      //  compute First Collided Method simulation time
      double first_collided_estimated_time = omp_get_wtime();
      fmt::print(
        "  ======   ==========   ================================\n\n");
      fmt::print("First Collided Method run time = {:.4f} \n\n",
        first_collided_estimated_time - end_volume_estimation);
    }
  }

  // Random ray power iteration loop
  while (simulation::current_batch < settings::n_batches) {

    // Initialize the current batch
    initialize_batch();
    initialize_generation();

    // Reset total starting particle weight used for normalizing tallies
    simulation::total_weight = 1.0;

    // Update external source if FIRST_COLLIDED_METHOD is used

    if (domain_->first_collided_mode) {
      if (settings::volume_online_option == 2) {
        domain_->update_external_flat_source();
      } else if (settings::volume_online_option == 3) {
        domain_->update_external_flat_source();
        domain_->update_external_linear_source();
      }
    }
    // Update source term (scattering + fission)
    domain_->update_neutron_source(k_eff_);

    // Reset scalar fluxes, iteration volume tallies, and region hit flags to
    // zero
    domain_->batch_reset();

    // Start timer for transport
    simulation::time_transport.start();

// Transport sweep over all random rays for the iteration
#pragma omp parallel for schedule(dynamic)                                     \
  reduction(+ : total_geometric_intersections_)
    for (int i = 0; i < simulation::work_per_rank; i++) {
      RandomRay ray(i, domain_.get());
      total_geometric_intersections_ +=
        ray.transport_history_based_single_ray();
    }

    simulation::time_transport.stop();

    // If using multiple MPI ranks, perform all reduce on all transport results
    domain_->all_reduce_replicated_source_regions();

    // Normalize scalar flux and update volumes
    domain_->normalize_scalar_flux_and_volumes(
      settings::n_particles * RandomRay::distance_active_);

    // Add source to scalar flux, compute number of FSR hits
    int64_t n_hits = domain_->add_source_to_scalar_flux();

    if (settings::run_mode == RunMode::EIGENVALUE) {
      // Compute random ray k-eff
      k_eff_ = domain_->compute_k_eff(k_eff_);

      // Store random ray k-eff into OpenMC's native k-eff variable
      global_tally_tracklength = k_eff_;
    }

    // Execute all tallying tasks, if this is an active batch
    if (simulation::current_batch > settings::n_inactive && mpi::master) {

      // Generate mapping between source regions and tallies
      if (!domain_->mapped_all_tallies_) {
        domain_->convert_source_regions_to_tallies();
      }

      // Use above mapping to contribute FSR flux data to appropriate tallies
      domain_->random_ray_tally();

      // Add this iteration's scalar flux estimate to final accumulated estimate
      domain_->accumulate_iteration_flux();
    }

    // Set phi_old = phi_new
    domain_->flux_swap();

    // Check for any obvious insabilities/nans/infs
    instability_check(n_hits, k_eff_, avg_miss_rate_);

    // Finalize the current batch
    finalize_generation();
    finalize_batch();
  } // End random ray power iteration loop

  if (domain_->first_collided_mode) {
    if (settings::volume_online_option != 1) {
      domain_->update_volume_uncollided_flux();
    }
  }
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

void RandomRaySimulation::output_simulation_results() const
{
  // Print random ray results
  if (mpi::master) {
    print_results_random_ray(total_geometric_intersections_,
      avg_miss_rate_ / settings::n_batches, negroups_,
      domain_->n_source_regions_, domain_->n_external_source_regions_);
    if (model::plots.size() > 0) {
      domain_->output_to_vtk();
    }
  }
}

// Apply a few sanity checks to catch obvious cases of numerical instability.
// Instability typically only occurs if ray density is extremely low.
void RandomRaySimulation::instability_check(
  int64_t n_hits, double k_eff, double& avg_miss_rate) const
{
  double percent_missed = ((domain_->n_source_regions_ - n_hits) /
                            static_cast<double>(domain_->n_source_regions_)) *
                          100.0;
  avg_miss_rate += percent_missed;

  if (mpi::master) {
    if (percent_missed > 10.0) {
      warning(fmt::format(
        "Very high FSR miss rate detected ({:.3f}%). Instability may occur. "
        "Increase ray density by adding more rays and/or active distance.",
        percent_missed));
    } else if (percent_missed > 0.01) {
      warning(
        fmt::format("Elevated FSR miss rate detected ({:.3f}%). Increasing "
                    "ray density by adding more rays and/or active "
                    "distance may improve simulation efficiency.",
          percent_missed));
    }

    if (k_eff > 10.0 || k_eff < 0.01 || !(std::isfinite(k_eff))) {
      fatal_error("Instability detected");
    }
  }
}

// Print random ray simulation results
void RandomRaySimulation::print_results_random_ray(
  uint64_t total_geometric_intersections, double avg_miss_rate, int negroups,
  int64_t n_source_regions, int64_t n_external_source_regions) const
{
  using namespace simulation;

  if (settings::verbosity >= 6) {
    double total_integrations = total_geometric_intersections * negroups;
    double time_per_integration =
      simulation::time_transport.elapsed() / total_integrations;
    double misc_time = time_total.elapsed() - time_update_src.elapsed() -
                       time_transport.elapsed() - time_tallies.elapsed() -
                       time_bank_sendrecv.elapsed();

    header("Simulation Statistics", 4);
    fmt::print(
      " Total Iterations                  = {}\n", settings::n_batches);
    fmt::print(" Flat Source Regions (FSRs)        = {}\n", n_source_regions);
    if (!domain_->first_collided_mode) {
      fmt::print(
        " FSRs Containing External Sources  = {}\n", n_external_source_regions);
    }
    fmt::print(" Total Geometric Intersections     = {:.4e}\n",
      static_cast<double>(total_geometric_intersections));
    fmt::print("   Avg per Iteration               = {:.4e}\n",
      static_cast<double>(total_geometric_intersections) / settings::n_batches);
    fmt::print("   Avg per Iteration per FSR       = {:.2f}\n",
      static_cast<double>(total_geometric_intersections) /
        static_cast<double>(settings::n_batches) / n_source_regions);
    fmt::print(" Avg FSR Miss Rate per Iteration   = {:.4f}%\n", avg_miss_rate);
    fmt::print(" Energy Groups                     = {}\n", negroups);
    fmt::print(
      " Total Integrations                = {:.4e}\n", total_integrations);
    fmt::print("   Avg per Iteration               = {:.4e}\n",
      total_integrations / settings::n_batches);

    header("Timing Statistics", 4);
    show_time("Total time for initialization", time_initialize.elapsed());
    show_time("Reading cross sections", time_read_xs.elapsed(), 1);
    show_time("Total simulation time", time_total.elapsed());
    show_time("Transport sweep only", time_transport.elapsed(), 1);
    show_time("Source update only", time_update_src.elapsed(), 1);
    show_time("Tally conversion only", time_tallies.elapsed(), 1);
    show_time("MPI source reductions only", time_bank_sendrecv.elapsed(), 1);
    show_time("Other iteration routines", misc_time, 1);
    if (settings::run_mode == RunMode::EIGENVALUE) {
      show_time("Time in inactive batches", time_inactive.elapsed());
    }
    show_time("Time in active batches", time_active.elapsed());
    show_time("Time writing statepoints", time_statepoint.elapsed());
    show_time("Total time for finalization", time_finalize.elapsed());
    show_time("Time per integration", time_per_integration);
  }

  if (settings::verbosity >= 4 && settings::run_mode == RunMode::EIGENVALUE) {
    header("Results", 4);
    fmt::print(" k-effective                       = {:.5f} +/- {:.5f}\n",
      simulation::keff, simulation::keff_std);
  }
}

} // namespace openmc
