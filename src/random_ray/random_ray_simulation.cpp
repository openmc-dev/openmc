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
#include "openmc/weight_windows.h"

namespace openmc {

//==============================================================================
// Non-member functions
//==============================================================================

void openmc_run_random_ray()
{
  //////////////////////////////////////////////////////////
  // Run forward simulation
  //////////////////////////////////////////////////////////

  // Check if adjoint calculation is needed. If it is, we will run the forward
  // calculation first and then the adjoint calculation later.
  bool adjoint_needed = FlatSourceDomain::adjoint_;

  // TODO: no initial condition needed fof td fixed source simulations
  // Check if this simulation is to establish an initial condition
  if (settings::kinetic_simulation) {
    simulation::is_initial_condition = true;
  }

  // Configure the domain for forward simulation
  FlatSourceDomain::adjoint_ = false;

  // If we're going to do a kinetic simulation, report that this is
  // the initial condition.
  if (settings::kinetic_simulation && mpi::master)
    header("KINETIC SIMULATION INITIAL CONDITION", 3);

  // If we're going to do an adjoint simulation afterwards, report that this is
  // the initial forward flux solve.
  if (adjoint_needed && mpi::master)
    header("FORWARD FLUX SOLVE", 3);

  // Initialize OpenMC general data structures
  openmc_simulation_init();

  // Validate that inputs meet requirements for random ray mode
  if (mpi::master)
    validate_random_ray_inputs();

  // Initialize Random Ray Simulation Object
  RandomRaySimulation sim;

  // Initialize fixed sources, if present
  sim.apply_fixed_sources_and_mesh_domains();

  // Begin main simulation timer
  simulation::time_total.start();

  // Execute random ray simulation
  sim.simulate();

  // End main simulation timer
  simulation::time_total.stop();

  // Normalize and save the final forward quantities
  sim.domain()->normalize_final_quantities();

  // Finalize OpenMC
  openmc_simulation_finalize();

  // Output all simulation results
  sim.output_simulation_results();

  if (settings::kinetic_simulation) {
    // Now do a second steady state simulation to correct the batch wise k-eff
    // estimates
    simulation::k_eff_correction = true;

    if (settings::kinetic_simulation && mpi::master)
      header("KINETIC SIMULATION INITIAL CONDITION (K-EFF CORRECTION)", 3);

    // If we're going to do an adjoint simulation afterwards, report that this
    // is the initial forward flux solve.
    if (adjoint_needed && mpi::master)
      header("FORWARD FLUX SOLVE", 3);

    // Initialize OpenMC general data structures
    openmc_simulation_init();

    double initial_k_eff = simulation::keff;

    sim.domain()->k_eff_ = initial_k_eff;
    sim.domain()->source_regions_.adjoint_reset();
    sim.domain()->propagate_final_quantities();
    sim.domain()->source_regions_.time_step_reset();

    // Begin main simulation timer
    simulation::time_total.start();

    // Execute random ray simulation
    sim.simulate();

    // End main simulation timer
    simulation::time_total.stop();

    // Normalize and save the final forward quantities
    sim.domain()->normalize_final_quantities();

    // Finalize OpenMC
    openmc_simulation_finalize();

    // Output all simulation results
    sim.output_simulation_results();

    //-------------------------------------------------------------------------
    // KINETIC SIMULATION

    warning(
      "Time-dependent explicit void treatment has not yet been "
      "implemented. Use caution when interpreting results from models with "
      "voids, as they may contain large inaccuracies.");
    sim.domain()->store_time_step_quantities(false);
    rename_statepoint_file(0);
    if (settings::output_tallies) {
      rename_tallies_file(0);
    }
    set_time_dependent_settings();

    // Timestepping loop
    // TODO: Add support for time-dependent restart
    for (int i = 0; i < settings::n_timesteps; i++) {
      simulation::current_timestep = i + 1;

      // Print simulation information
      if (mpi::master) {
        std::string message = fmt::format(
          "KINETIC SIMULATION TIME STEP {0}", simulation::current_timestep);
        const char* msg = message.c_str();
        header(msg, 3);
      }

      // If we're going to do an adjoint simulation afterwards, report that this
      // is the initial forward flux solve.
      if (adjoint_needed && mpi::master)
        header("FORWARD FLUX SOLVE", 3);

      reset_timers();

      // Initialize OpenMC general data structures
      openmc_simulation_init();

      sim.domain()->k_eff_ = initial_k_eff;
      sim.domain()->source_regions_.adjoint_reset();
      sim.domain()->propagate_final_quantities();
      sim.domain()->source_regions_.time_step_reset();

      // Compute RHS backward differences to be used later
      sim.domain()->compute_rhs_bd_quantities();

      // Update time dependent cross section based on the density
      sim.domain()->update_material_density(i);

      // Begin main simulation timer
      simulation::time_total.start();

      // Execute random ray simulation
      sim.simulate();

      // End main simulation timer
      simulation::time_total.stop();

      // Finalize OpenMC
      openmc_simulation_finalize();

      // Output all simulation results
      sim.output_simulation_results();

      // Rename statepoint and tallies file
      rename_statepoint_file(simulation::current_timestep);
      if (settings::output_tallies) {
        rename_tallies_file(simulation::current_timestep);
      }

      // Normalize and store final quantities for next time step
      sim.domain()->normalize_final_quantities();
      sim.domain()->store_time_step_quantities();

      // Advance time
      simulation::current_time += settings::dt;
    }
  }

  //////////////////////////////////////////////////////////
  // Run adjoint simulation (if enabled)
  //////////////////////////////////////////////////////////

  if (!adjoint_needed) {
    return;
  }

  reset_timers();

  // Configure the domain for adjoint simulation
  FlatSourceDomain::adjoint_ = true;

  if (mpi::master)
    header("ADJOINT FLUX SOLVE", 3);

  sim.domain()->k_eff_ = 1.0;

  // Initialize adjoint fixed sources, if present
  sim.prepare_fixed_sources_adjoint();

  // Transpose scattering matrix
  sim.domain()->transpose_scattering_matrix();

  // Swap nu_sigma_f and chi
  sim.domain()->nu_sigma_f_.swap(sim.domain()->chi_);

  // Initialize OpenMC general data structures
  openmc_simulation_init();

  // Begin main simulation timer
  simulation::time_total.start();

  // Execute random ray simulation
  sim.simulate();

  // End main simulation timer
  simulation::time_total.stop();

  // Finalize OpenMC
  openmc_simulation_finalize();

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

      // TODO: add support for prompt and delayed fission
      case SCORE_PRECURSORS: {
        if (settings::kinetic_simulation && settings::create_delayed_neutrons) {
          break;
        } else {
          fatal_error("Invalid score specified in tallies.xml. Precursors can "
                      "only be scored for kinetic simulations.");
        }
      }
      default:
        fatal_error(
          "Invalid score specified. Only flux, total, fission, nu-fission, and "
          "event scores are supported in random ray mode. (precursors are "
          "supported for "
          "kinetic simulations when delayed neutrons are turned on).");
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
      case FilterType::PARTICLE:
        break;
      case FilterType::DELAYED_GROUP:
        if (settings::kinetic_simulation) {
          break;
        } else {
          fatal_error("Invalid filter specified in tallies.xml. Kinetic "
                      "simulations is required "
                      "to tally with a delayed_group filter.");
        }
      default:
        fatal_error("Invalid filter specified. Only cell, cell_instance, "
                    "distribcell, energy, material, mesh, and universe filters "
                    "are supported in random ray mode (delayed_group is "
                    "supported for kinetic simulations).");
      }
    }
  }

  // TODO: validate kinetic data is present
  //  Validate MGXS data
  ///////////////////////////////////////////////////////////////////
  for (auto& material : data::mg.macro_xs_) {
    if (!material.is_isotropic) {
      fatal_error("Anisotropic MGXS detected. Only isotropic XS data sets "
                  "supported in random ray mode.");
    }
    if (material.get_xsdata().size() > 1) {
      warning("Non-isothermal MGXS detected. Only isothermal XS data sets "
              "supported in random ray mode. Using lowest temperature.");
    }
    for (int g = 0; g < data::mg.num_energy_groups_; g++) {
      if (material.exists_in_model) {
        // Temperature and angle indices, if using multiple temperature
        // data sets and/or anisotropic data sets.
        // TODO: Currently assumes we are only using single temp/single angle
        // data.
        const int t = 0;
        const int a = 0;
        double sigma_t =
          material.get_xs(MgxsType::TOTAL, g, NULL, NULL, NULL, t, a);
        if (sigma_t <= 0.0) {
          fatal_error("No zero or negative total macroscopic cross sections "
                      "allowed in random ray mode. If the intention is to make "
                      "a void material, use a cell fill of 'None' instead.");
        }
      }
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

      // Validate that a domain ID was specified OR that it is a point source
      auto sp = dynamic_cast<SpatialPoint*>(is->space());
      if (is->domain_ids().size() == 0 && !sp) {
        fatal_error("Fixed sources must be point source or spatially "
                    "constrained by domain id (cell, material, or universe) in "
                    "random ray mode.");
      } else if (is->domain_ids().size() > 0 && sp) {
        // If both a domain constraint and a non-default point source location
        // are specified, notify user that domain constraint takes precedence.
        if (sp->r().x == 0.0 && sp->r().y == 0.0 && sp->r().z == 0.0) {
          warning("Fixed source has both a domain constraint and a point "
                  "type spatial distribution. The domain constraint takes "
                  "precedence in random ray mode -- point source coordinate "
                  "will be ignored.");
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
    const auto& openmc_plottable = model::plots[p];
    Plot* openmc_plot = dynamic_cast<Plot*>(openmc_plottable.get());

    // Random ray plots only support voxel plots
    if (!openmc_plot) {
      warning(fmt::format(
        "Plot {} will not be used for end of simulation data plotting -- only "
        "voxel plotting is allowed in random ray mode.",
        openmc_plottable->id()));
      continue;
    } else if (openmc_plot->type_ != Plot::PlotType::voxel) {
      warning(fmt::format(
        "Plot {} will not be used for end of simulation data plotting -- only "
        "voxel plotting is allowed in random ray mode.",
        openmc_plottable->id()));
      continue;
    }
  }

  // Warn about slow MPI domain replication, if detected
  ///////////////////////////////////////////////////////////////////
#ifdef OPENMC_MPI
  if (mpi::n_procs > 1) {
    warning(
      "MPI parallelism is not supported by the random ray solver. All work "
      "will be performed by rank 0. Domain decomposition may be implemented in "
      "the future to provide efficient MPI scaling.");
  }
#endif

  // Warn about instability resulting from linear sources in small regions
  // when generating weight windows with FW-CADIS and an overlaid mesh.
  ///////////////////////////////////////////////////////////////////
  if (RandomRay::source_shape_ == RandomRaySourceShape::LINEAR &&
      variance_reduction::weight_windows.size() > 0) {
    warning(
      "Linear sources may result in negative fluxes in small source regions "
      "generated by mesh subdivision. Negative sources may result in low "
      "quality FW-CADIS weight windows. We recommend you use flat source mode "
      "when generating weight windows with an overlaid mesh tally.");
  }
}

void openmc_reset_random_ray()
{
  FlatSourceDomain::volume_estimator_ = RandomRayVolumeEstimator::HYBRID;
  FlatSourceDomain::volume_normalized_flux_tallies_ = false;
  FlatSourceDomain::adjoint_ = false;
  FlatSourceDomain::mesh_domain_map_.clear();
  RandomRay::ray_source_.reset();
  RandomRay::source_shape_ = RandomRaySourceShape::FLAT;
  RandomRay::sample_method_ = RandomRaySampleMethod::PRNG;
}

void write_random_ray_hdf5(hid_t group)
{
  hid_t random_ray_group = create_group(group, "random_ray");
  switch (RandomRay::source_shape_) {
  case RandomRaySourceShape::FLAT:
    write_dataset(random_ray_group, "source_shape", "flat");
    break;
  case RandomRaySourceShape::LINEAR:
    write_dataset(random_ray_group, "source_shape", "linear");
    break;
  case RandomRaySourceShape::LINEAR_XY:
    write_dataset(random_ray_group, "source_shape", "linear xy");
    break;
  default:
    break;
  }

  switch (FlatSourceDomain::volume_estimator_) {
  case RandomRayVolumeEstimator::SIMULATION_AVERAGED:
    write_dataset(random_ray_group, "volume_estimator", "simulation averaged");
    break;
  case RandomRayVolumeEstimator::NAIVE:
    write_dataset(random_ray_group, "volume_estimator", "naive");
    break;
  case RandomRayVolumeEstimator::HYBRID:
    write_dataset(random_ray_group, "volume_estimator", "hybrid");
    break;
  default:
    break;
  }

  write_dataset(
    random_ray_group, "distance_active", RandomRay::distance_active_);
  write_dataset(
    random_ray_group, "distance_inactive", RandomRay::distance_inactive_);
  write_dataset(random_ray_group, "volume_normalized_flux_tallies",
    FlatSourceDomain::volume_normalized_flux_tallies_);
  write_dataset(random_ray_group, "adjoint_mode", FlatSourceDomain::adjoint_);

  write_dataset(random_ray_group, "avg_miss_rate", RandomRay::avg_miss_rate_);
  write_dataset(
    random_ray_group, "n_source_regions", RandomRay::n_source_regions_);
  write_dataset(random_ray_group, "n_external_source_regions",
    RandomRay::n_external_source_regions_);
  write_dataset(random_ray_group, "n_geometric_intersections",
    RandomRay::total_geometric_intersections_);
  int64_t n_integrations =
    RandomRay::total_geometric_intersections_ * data::mg.num_energy_groups_;
  write_dataset(random_ray_group, "n_integrations", n_integrations);

  if (settings::kinetic_simulation && !simulation::is_initial_condition) {
    write_dataset(random_ray_group, "bd_order", RandomRay::bd_order_);
    switch (RandomRay::time_method_) {
    case RandomRayTimeMethod::ISOTROPIC:
      write_dataset(random_ray_group, "time_method", "isotropic");
      break;
    case RandomRayTimeMethod::PROPAGATION:
      write_dataset(random_ray_group, "time_method", "propogation");
      break;
    default:
      break;
    }
  }
  close_group(random_ray_group);
}

//-----------------------------------------------------------------------------
// Non-member functions for kinetic simulations

void set_time_dependent_settings()
{
  // Reset flags
  simulation::is_initial_condition = false;
  simulation::k_eff_correction = false;

  // Set current time
  simulation::current_time = settings::dt;
}

// TODO: condense this into one function with rename_tallies_file and use char
// or string arguments
void rename_statepoint_file(int i)
{
  // Rename statepoint file
  std::string old_filename_ = fmt::format(
    "{0}statepoint.{1}.h5", settings::path_output, settings::n_batches);
  std::string new_filename_ =
    fmt::format("{0}openmc_td_simulation_{1}.h5", settings::path_output, i);

  const char* old_fname = old_filename_.c_str();
  const char* new_fname = new_filename_.c_str();
  std::rename(old_fname, new_fname);
}

void rename_tallies_file(int i)
{
  // Rename tallies file
  std::string old_filename_ =
    fmt::format("{0}tallies.out", settings::path_output);
  std::string new_filename_ =
    fmt::format("{0}tallies_{1}.out", settings::path_output, i);

  const char* old_fname = old_filename_.c_str();
  const char* new_fname = new_filename_.c_str();
  std::rename(old_fname, new_fname);
}

//==============================================================================
// RandomRaySimulation implementation
//==============================================================================

RandomRaySimulation::RandomRaySimulation()
  : negroups_(data::mg.num_energy_groups_),
    ndgroups_(data::mg.num_delayed_groups_)
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
  case RandomRaySourceShape::LINEAR_XY:
    domain_ = make_unique<LinearSourceDomain>();
    break;
  default:
    fatal_error("Unknown random ray source shape");
  }

  // Convert OpenMC native MGXS into a more efficient format
  // internal to the random ray solver
  domain_->flatten_xs();

  // Initialize vectors used for the steady state
  if (settings::kinetic_simulation) {
    static_k_eff_;
    static_fission_rate_;
  }
}

void RandomRaySimulation::apply_fixed_sources_and_mesh_domains()
{
  domain_->apply_meshes();
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    // Transfer external source user inputs onto random ray source regions
    domain_->convert_external_sources();
    domain_->count_external_source_regions();
  }
}

void RandomRaySimulation::prepare_fixed_sources_adjoint()
{
  domain_->source_regions_.adjoint_reset();
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    domain_->set_adjoint_sources();
  }
}

void RandomRaySimulation::simulate()
{
  // Random ray power iteration loop
  while (simulation::current_batch < settings::n_batches) {
    // Initialize the current batch
    initialize_batch();
    initialize_generation();

    // MPI not supported in random ray solver, so all work is done by rank 0
    // TODO: Implement domain decomposition for MPI parallelism
    if (mpi::master) {

      // Reset total starting particle weight used for normalizing tallies
      simulation::total_weight = 1.0;

      // TODO: add update source convenience function
      // domain_->compute_neutron_source()
      // Update source term (scattering + fission)
      domain_->update_all_neutron_sources();

      // Reset scalar fluxes, iteration volume tallies, and region hit flags
      // to zero
      domain_->batch_reset();

      // At the beginning of the simulation, if mesh subdivision is in use, we
      // need to swap the main source region container into the base container,
      // as the main source region container will be used to hold the true
      // subdivided source regions. The base container will therefore only
      // contain the external source region information, the mesh indices,
      // material properties, and initial guess values for the flux/source.

      // Start timer for transport
      simulation::time_transport.start();

// Transport sweep over all random rays for the iteration
#pragma omp parallel for schedule(dynamic)                                     \
  reduction(+ : total_geometric_intersections_)
      for (int i = 0; i < settings::n_particles; i++) {
        RandomRay ray(i, domain_.get());
        total_geometric_intersections_ +=
          ray.transport_history_based_single_ray();
      }

      simulation::time_transport.stop();

      // Add any newly discovered source regions to the main source region
      // container.
      domain_->finalize_discovered_source_regions();

      // Normalize scalar flux and update volumes
      domain_->normalize_scalar_flux_and_volumes(
        settings::n_particles * RandomRay::distance_active_);

      // Add source to scalar flux, compute number of FSR hits
      int64_t n_hits = domain_->add_source_to_scalar_flux();

      // Apply transport stabilization factors
      domain_->apply_transport_stabilization();

      if (settings::run_mode == RunMode::EIGENVALUE) {
        // Compute random ray k-eff
        if (!settings::kinetic_simulation ||
            settings::kinetic_simulation && simulation::is_initial_condition) {
          domain_->compute_k_eff();
          if (simulation::k_eff_correction) {
            static_fission_rate_.push_back(domain_->fission_rate_);
            static_k_eff_.push_back(domain_->k_eff_);
          }
        } else {
          domain_->k_eff_ = static_k_eff_[simulation::current_batch - 1];
          domain_->fission_rate_ =
            static_fission_rate_[simulation::current_batch - 1];
        }

        // Store random ray k-eff into OpenMC's native k-eff variable
        global_tally_tracklength = domain_->k_eff_;
      }

      // Compute precursors if delayed neutrons are turned on
      if (settings::kinetic_simulation && settings::create_delayed_neutrons)
        domain_->compute_all_precursors();

      // Execute all tallying tasks, if this is an active batch
      if (simulation::current_batch > settings::n_inactive) {

        // Add this iteration's scalar flux estimate to final accumulated
        // estimate
        domain_->accumulate_iteration_quantities();

        // Use above mapping to contribute FSR flux data to appropriate
        // tallies
        domain_->random_ray_tally();
      }

      // Set phi_old = phi_new
      domain_->flux_swap();
      if (settings::kinetic_simulation && settings::create_delayed_neutrons) {
        domain_->precursors_swap();
      }

      // Check for any obvious insabilities/nans/infs
      instability_check(n_hits, domain_->k_eff_, avg_miss_rate_);
    } // End MPI master work

    // Store simulation metrics
    RandomRay::avg_miss_rate_ = avg_miss_rate_ / settings::n_batches;
    RandomRay::total_geometric_intersections_ = total_geometric_intersections_;
    RandomRay::n_external_source_regions_ = domain_->n_external_source_regions_;
    RandomRay::n_source_regions_ = domain_->n_source_regions();

    // Finalize the current batch
    finalize_generation();
    finalize_batch();
  } // End random ray power iteration loop

  domain_->count_external_source_regions();
}

void RandomRaySimulation::output_simulation_results() const
{
  // Print random ray results
  if (mpi::master) {
    print_results_random_ray();
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
  double percent_missed = ((domain_->n_source_regions() - n_hits) /
                            static_cast<double>(domain_->n_source_regions())) *
                          100.0;
  avg_miss_rate += percent_missed;

  if (mpi::master) {
    if (percent_missed > 10.0) {
      warning(fmt::format(
        "Very high FSR miss rate detected ({:.3f}%). Instability may occur. "
        "Increase ray density by adding more rays and/or active distance.",
        percent_missed));
    } else if (percent_missed > 1.0) {
      warning(
        fmt::format("Elevated FSR miss rate detected ({:.3f}%). Increasing "
                    "ray density by adding more rays and/or active "
                    "distance may improve simulation efficiency.",
          percent_missed));
    }

    if (k_eff > 10.0 || k_eff < 0.01 || !(std::isfinite(k_eff))) {
      fatal_error(fmt::format("Instability detected: k-eff = {:.5f}", k_eff));
    }
  }
}

// Print random ray simulation results
void RandomRaySimulation::print_results_random_ray() const
{
  using namespace simulation;

  if (settings::verbosity >= 6) {
    double total_integrations =
      RandomRay::total_geometric_intersections_ * negroups_;
    double time_per_integration =
      simulation::time_transport.elapsed() / total_integrations;
    double misc_time = time_total.elapsed() - time_update_src.elapsed() -
                       time_transport.elapsed() - time_tallies.elapsed() -
                       time_bank_sendrecv.elapsed();

    if (settings::kinetic_simulation && !simulation::is_initial_condition) {
      misc_time -= time_update_bd_vectors.elapsed();
    }
    header("Simulation Statistics", 4);
    fmt::print(
      " Total Iterations                  = {}\n", settings::n_batches);
    fmt::print(
      " Number of Rays per Iteration      = {}\n", settings::n_particles);
    fmt::print(" Inactive Distance                 = {} cm\n",
      RandomRay::distance_inactive_);
    fmt::print(" Active Distance                   = {} cm\n",
      RandomRay::distance_active_);
    fmt::print(" Source Regions (SRs)              = {}\n",
      RandomRay::n_source_regions_);
    fmt::print(" SRs Containing External Sources   = {}\n",
      RandomRay::n_external_source_regions_);
    fmt::print(" Total Geometric Intersections     = {:.4e}\n",
      static_cast<double>(RandomRay::total_geometric_intersections_));
    fmt::print("   Avg per Iteration               = {:.4e}\n",
      static_cast<double>(RandomRay::total_geometric_intersections_) /
        settings::n_batches);
    fmt::print("   Avg per Iteration per SR        = {:.2f}\n",
      static_cast<double>(RandomRay::total_geometric_intersections_) /
        static_cast<double>(settings::n_batches) /
        RandomRay::n_source_regions_);
    fmt::print(" Avg SR Miss Rate per Iteration    = {:.4f}%\n",
      RandomRay::avg_miss_rate_);
    fmt::print(" Energy Groups                     = {}\n", negroups_);
    if (settings::kinetic_simulation)
      fmt::print(" Delay Groups                      = {}\n", ndgroups_);
    fmt::print(
      " Total Integrations                = {:.4e}\n", total_integrations);
    fmt::print("   Avg per Iteration               = {:.4e}\n",
      total_integrations / settings::n_batches);

    std::string estimator;
    switch (domain_->volume_estimator_) {
    case RandomRayVolumeEstimator::SIMULATION_AVERAGED:
      estimator = "Simulation Averaged";
      break;
    case RandomRayVolumeEstimator::NAIVE:
      estimator = "Naive";
      break;
    case RandomRayVolumeEstimator::HYBRID:
      estimator = "Hybrid";
      break;
    default:
      fatal_error("Invalid volume estimator type");
    }
    fmt::print(" Volume Estimator Type             = {}\n", estimator);

    std::string adjoint_true = (FlatSourceDomain::adjoint_) ? "ON" : "OFF";
    fmt::print(" Adjoint Flux Mode                 = {}\n", adjoint_true);

    std::string shape;
    switch (RandomRay::source_shape_) {
    case RandomRaySourceShape::FLAT:
      shape = "Flat";
      break;
    case RandomRaySourceShape::LINEAR:
      shape = "Linear";
      break;
    case RandomRaySourceShape::LINEAR_XY:
      shape = "Linear XY";
      break;
    default:
      fatal_error("Invalid random ray source shape");
    }
    fmt::print(" Source Shape                      = {}\n", shape);
    std::string sample_method =
      (RandomRay::sample_method_ == RandomRaySampleMethod::PRNG) ? "PRNG"
                                                                 : "Halton";
    fmt::print(" Sample Method                     = {}\n", sample_method);

    if (domain_->is_transport_stabilization_needed_) {
      fmt::print(" Transport XS Stabilization Used   = YES (rho = {:.3f})\n",
        FlatSourceDomain::diagonal_stabilization_rho_);
    } else {
      fmt::print(" Transport XS Stabilization Used   = NO\n");
    }
    if (settings::kinetic_simulation && !simulation::is_initial_condition) {
      std::string time_method =
        (RandomRay::time_method_ == RandomRayTimeMethod::ISOTROPIC)
          ? "ISOTROPIC"
          : "PROPAGATION";
      fmt::print(" Time Method                       = {}\n", time_method);
      fmt::print(
        " Backwards Difference Order        = {}\n", RandomRay::bd_order_);
    }

    header("Timing Statistics", 4);
    show_time("Total time for initialization", time_initialize.elapsed());
    show_time("Reading cross sections", time_read_xs.elapsed(), 1);
    show_time("Total simulation time", time_total.elapsed());
    show_time("Transport sweep only", time_transport.elapsed(), 1);
    show_time("Source update only", time_update_src.elapsed(), 1);
    if (settings::kinetic_simulation && settings::create_delayed_neutrons) {
      show_time(
        "Precursor computation only", time_compute_precursors.elapsed(), 1);
      misc_time -= time_compute_precursors.elapsed();
    }
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
