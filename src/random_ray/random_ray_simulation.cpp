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
#include "openmc/random_ray/decomposition_map.h"
#include "openmc/random_ray/ray_bank.h"

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

  // Configure the domain for forward simulation
  FlatSourceDomain::adjoint_ = false;

  // If we're going to do an adjoint simulation afterwards, report that this is
  // the initial forward flux solve.
  if (adjoint_needed && mpi::master)
    header("FORWARD FLUX SOLVE", 3);

  // Initialize OpenMC general data structures
  openmc_simulation_init();

  // Validate that inputs meet requirements for random ray mode
  if (mpi::master)
    validate_random_ray_inputs();

  // Declare forward flux so that it can be saved for later adjoint simulation
  vector<double> forward_flux;
  SourceRegionContainer forward_source_regions;
  SourceRegionContainer forward_base_source_regions;
  std::unordered_map<SourceRegionKey, int64_t, SourceRegionKey::HashFunctor>
    forward_source_region_map;

  {
    // Initialize Random Ray Simulation Object
    RandomRaySimulation sim;

    // Initialize fixed sources, if present
    sim.apply_fixed_sources_and_mesh_domains();

    // Initialize subdomains for MPI ranks
    simulation::time_decomposition_handling.start();
    mpi::decomp_map.generate_rank_centers();
    // mpi::decomp_map.create(sim.domain());
    simulation::time_decomposition_handling.stop();

    // Begin main simulation timer
    simulation::time_total.start();

    // Execute random ray simulation
    sim.simulate();

    // End main simulation timer
    simulation::time_total.stop();

    // Normalize and save the final forward flux
    sim.domain()->serialize_final_fluxes(forward_flux);

    double source_normalization_factor =
      sim.domain()->compute_fixed_source_normalization_factor() /
      (settings::n_batches - settings::n_inactive);

#pragma omp parallel for
    for (uint64_t i = 0; i < forward_flux.size(); i++) {
      forward_flux[i] *= source_normalization_factor;
    }

    forward_source_regions = sim.domain()->source_regions_;
    forward_source_region_map = sim.domain()->source_region_map_;
    forward_base_source_regions = sim.domain()->base_source_regions_;

    // Finalize OpenMC
    openmc_simulation_finalize();

    // Output all simulation results
    sim.output_simulation_results();
  }

  //////////////////////////////////////////////////////////
  // Run adjoint simulation (if enabled)
  //////////////////////////////////////////////////////////

  if (adjoint_needed) {
    reset_timers();

    // Configure the domain for adjoint simulation
    FlatSourceDomain::adjoint_ = true;

    if (mpi::master)
      header("ADJOINT FLUX SOLVE", 3);

    // Initialize OpenMC general data structures
    openmc_simulation_init();

    // Initialize Random Ray Simulation Object
    RandomRaySimulation adjoint_sim;

    // Initialize adjoint fixed sources, if present
    adjoint_sim.prepare_fixed_sources_adjoint(forward_flux,
      forward_source_regions, forward_base_source_regions,
      forward_source_region_map);

    // Transpose scattering matrix
    adjoint_sim.domain()->transpose_scattering_matrix();

    // Swap nu_sigma_f and chi
    adjoint_sim.domain()->nu_sigma_f_.swap(adjoint_sim.domain()->chi_);

    // Begin main simulation timer
    simulation::time_total.start();

    // Execute random ray simulation
    adjoint_sim.simulate();

    // End main simulation timer
    simulation::time_total.stop();

    // Finalize OpenMC
    openmc_simulation_finalize();

    // Output all simulation results
    adjoint_sim.output_simulation_results();
  }

  printf("Finished adjoint simulation\n");

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
      case FilterType::PARTICLE:
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

  // Warn about instability resulting from linear sources in small regions
  // when generating weight windows with FW-CADIS and an overlaid mesh.
  ///////////////////////////////////////////////////////////////////
  if (RandomRay::source_shape_ == RandomRaySourceShape::LINEAR &&
      RandomRay::mesh_subdivision_enabled_ &&
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
  RandomRay::mesh_subdivision_enabled_ = false;
  RandomRay::sample_method_ = RandomRaySampleMethod::PRNG;
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
  case RandomRaySourceShape::LINEAR_XY:
    domain_ = make_unique<LinearSourceDomain>();
    break;
  default:
    fatal_error("Unknown random ray source shape");
  }

  // Convert OpenMC native MGXS into a more efficient format
  // internal to the random ray solver
  domain_->flatten_xs();
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

void RandomRaySimulation::prepare_fixed_sources_adjoint(
  vector<double>& forward_flux, SourceRegionContainer& forward_source_regions,
  SourceRegionContainer& forward_base_source_regions,
  std::unordered_map<SourceRegionKey, int64_t, SourceRegionKey::HashFunctor>&
    forward_source_region_map)
{
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    if (RandomRay::mesh_subdivision_enabled_) {
      domain_->source_regions_ = forward_source_regions;
      domain_->source_region_map_ = forward_source_region_map;
      domain_->base_source_regions_ = forward_base_source_regions;
      domain_->source_regions_.adjoint_reset();
    }
    domain_->set_adjoint_sources(forward_flux);
  }
}

void RandomRaySimulation::simulate()
{
  // Create ray bank
  RayBank RB;

  // Random ray power iteration loop
  while (simulation::current_batch < settings::n_batches) {
    // Initialize the current batch
    initialize_batch();
    initialize_generation();

    // Reset total starting particle weight used for normalizing tallies
    simulation::total_weight = 1.0;

    // Update source term (scattering + fission)
    domain_->update_neutron_source(k_eff_);

    // Reset scalar fluxes, iteration volume tallies, and region hit flags to
    // zero
    domain_->batch_reset();

    // At the beginning of the simulation, if mesh subvivision is in use, we
    // need to swap the main source region container into the base container,
    // as the main source region container will be used to hold the true
    // subdivided source regions. The base container will therefore only
    // contain the external source region information, the mesh indices,
    // material properties, and initial guess values for the flux/source.
    if (RandomRay::mesh_subdivision_enabled_ &&
        simulation::current_batch == 1 && !FlatSourceDomain::adjoint_) {
      domain_->prepare_base_source_regions();
    }

    // printf("Rank %d Test1 \n", mpi::rank);

    // Transport sweep over all random rays for the iteration
    if (mpi::n_procs > 1){

      // printf("1: Rank %d Starting transport sweep \n", mpi::rank);
      transport_sweep_decomp(RB);
      // printf("2: Rank %d Ended transport sweep \n", mpi::rank);

      // bool test  = domain_->discovered_source_regions_.contains(SourceRegionKey(4, 104171));
      // printf("2 RANK %d: Source region (4, 104171) discovered? %d \n", mpi::rank, test);

      // bool test2 = false;
      // if(domain_->source_region_map_.find(SourceRegionKey(4, 104171)) != domain_->source_region_map_.end()){
      //   test2 = true;
      // }
      // printf("2 RANK %d: Source region (4, 104171) in main map? %d \n", mpi::rank, test2);
      // if (test) {
      //   SourceRegion& sr = domain_->discovered_source_regions_[SourceRegionKey(4, 104171)];
      //   printf("RANK %d: volume: %f\n", mpi::rank, sr.scalars_.volume_);
      // }

      // Update decomposition map with newly discovered source regions
      simulation::time_decomposition_handling.start();
      mpi::decomp_map.update(domain_->discovered_source_regions_);
      simulation::time_decomposition_handling.stop();
      // printf("3: Rank %d Updated decomp map \n", mpi::rank);

    } else {
      transport_sweep();
    }

    // printf("Rank %d Test6 \n", mpi::rank);

    // bool test  = domain_->discovered_source_regions_.contains(SourceRegionKey(4, 104171));
    // printf("3 RANK %d: Source region (4, 104171) discovered? %d \n", mpi::rank, test);
    // if (test) {
    //   SourceRegion& sr = domain_->discovered_source_regions_[SourceRegionKey(4, 104171)];
    //   printf("RANK %d: volume: %f\n", mpi::rank, sr.scalars_.volume_);
    // }


    // bool test2 = false;
    // if(domain_->source_region_map_.find(SourceRegionKey(4, 104171)) != domain_->source_region_map_.end()){
    //   test2 = true;
    // }
    // printf("3 RANK %d: Source region (4, 104171) in main map? %d \n", mpi::rank, test2);

    // If using mesh subdivision, add any newly discovered source regions
    // to the main source region container.
    if (RandomRay::mesh_subdivision_enabled_) {
      domain_->finalize_discovered_source_regions();
    }

    // test  = domain_->discovered_source_regions_.contains(SourceRegionKey(4, 104171));
    // printf("4 RANK %d: Source region (4, 104171) discovered? %d \n", mpi::rank, test);

    // test2 = false;
    // if(domain_->source_region_map_.find(SourceRegionKey(4, 104171)) != domain_->source_region_map_.end()){
    //   test2 = true;
    // }
    // printf("4 RANK %d: Source region (4, 104171) in main map? %d \n", mpi::rank, test2);

    // Normalize scalar flux and update volumes
    domain_->normalize_scalar_flux_and_volumes(
      settings::n_particles * RandomRay::distance_active_);

    // Add source to scalar flux, compute number of FSR hits
    int64_t n_hits = domain_->add_source_to_scalar_flux();

    // Apply transport stabilization factors
    domain_->apply_transport_stabilization();

    if (settings::run_mode == RunMode::EIGENVALUE) {
      // Compute random ray k-eff
      k_eff_ = domain_->compute_k_eff(k_eff_);

      // Store random ray k-eff into OpenMC's native k-eff variable
      global_tally_tracklength = k_eff_;
    }

    // Execute all tallying tasks, if this is an active batch
    if (simulation::current_batch > settings::n_inactive) {

      // Add this iteration's scalar flux estimate to final accumulated
      // estimate
      domain_->accumulate_iteration_flux();

      // Generate mapping between source regions and tallies
      if (!domain_->mapped_all_tallies_ &&
          !RandomRay::mesh_subdivision_enabled_) {
        domain_->convert_source_regions_to_tallies(0);
      }

      // Use above mapping to contribute FSR flux data to appropriate
      // tallies
      domain_->random_ray_tally();
    }

    // Set phi_old = phi_new
    domain_->flux_swap();

    // Check for any obvious insabilities/nans/infs
    instability_check(n_hits, k_eff_, avg_miss_rate_);

    // Finalize the current batch
    finalize_generation();
    finalize_batch();
  } // End random ray power iteration loop

  if (mpi::n_procs > 1){

    // Compute total intersections and load balancing data
    simulation::time_decomposition_handling.start();

    // Exchange intersection data
    uint64_t total_geometric_intersections_sum = 0;
    MPI_Allreduce(&total_geometric_intersections_, &total_geometric_intersections_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi::intracomm);

    // Local work load 
    float rank_load_local = (total_geometric_intersections_/static_cast<double>(total_geometric_intersections_sum))*100.0;

    total_geometric_intersections_ = total_geometric_intersections_sum;
    rank_load_.resize(mpi::n_procs);
    // send load data to master rank (rank 0)
    if (mpi::master) {
      rank_load_[0] = rank_load_local;
      
      for (int i = 1; i < mpi::n_procs; i++) {
        MPI_Recv(&rank_load_[i], 1, MPI_FLOAT, i, 1, mpi::intracomm, MPI_STATUS_IGNORE);
      }
    } else {
      MPI_Send(&rank_load_local, 1, MPI_FLOAT, 0, 1, mpi::intracomm);
    }

    // Average number of ray communications between ranks per batch
    avg_num_comms_ = avg_num_comms_/settings::n_batches; //TODO: should this be double?

    simulation::time_decomposition_handling.stop();
  }

  domain_->count_external_source_regions();
}

void RandomRaySimulation::output_simulation_results() const
{
  // Print random ray results
  if (mpi::master) {
    print_results_random_ray(total_geometric_intersections_,
      avg_miss_rate_ / settings::n_batches, negroups_,
      domain_->n_source_regions(), domain_->n_external_source_regions_, avg_num_comms_);
  }

  if (model::plots.size() > 0) {
      if (mpi::n_procs > 1){
        domain_->output_to_vtk_decomp();
      } else {
        domain_->output_to_vtk();
      }
    }
  
}

// Apply a few sanity checks to catch obvious cases of numerical instability.
// Instability typically only occurs if ray density is extremely low.
void RandomRaySimulation::instability_check(
  int64_t n_hits, double k_eff, double& avg_miss_rate) const
{

  uint64_t n_source_regions = domain_->n_source_regions();

  if (mpi::n_procs > 1){
    // Reduce n_hits and n_source_regions on master rank to compute
    // global miss rate
    // printf("RANK %d: n_hits: %ld, n_source_regions: %ld \n", mpi::rank, n_hits, n_source_regions);

    simulation::time_decomposition_handling.stop();
    if (mpi::master) {
      MPI_Reduce(MPI_IN_PLACE, &n_hits, 1, MPI_LONG_LONG, MPI_SUM, 0, mpi::intracomm);
      MPI_Reduce(MPI_IN_PLACE, &n_source_regions, 1, MPI_LONG_LONG, MPI_SUM, 0, mpi::intracomm);
    } else {
      MPI_Reduce(&n_hits, nullptr, 1, MPI_LONG_LONG, MPI_SUM, 0, mpi::intracomm);
      MPI_Reduce(&n_source_regions, nullptr, 1, MPI_LONG_LONG, MPI_SUM, 0, mpi::intracomm);
    }
    simulation::time_decomposition_handling.stop();
  }

  if (mpi::master) {
    // printf("TOTAL: n_hits: %ld, n_source_regions: %ld \n", n_hits, n_source_regions);

    double percent_missed = ((n_source_regions - n_hits) /
                            static_cast<double>(n_source_regions)) *
                          100.0;
    avg_miss_rate += percent_missed;

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
      fatal_error("Instability detected");
    }
  }
}

// Print random ray simulation results
void RandomRaySimulation::print_results_random_ray(
  uint64_t total_geometric_intersections, double avg_miss_rate, int negroups,
  int64_t n_source_regions, int64_t n_external_source_regions, uint64_t avg_num_comms) const
{
  using namespace simulation;

  if (settings::verbosity >= 6) {
    double total_integrations = total_geometric_intersections * negroups;
    double time_per_integration =
      (time_transport.elapsed() - time_ray_buffering.elapsed())
      / total_integrations;
    double misc_time = time_total.elapsed() - time_update_src.elapsed() -
                       time_transport.elapsed() - time_tallies.elapsed() -
                       time_bank_sendrecv.elapsed();

    header("Simulation Statistics", 4);
    fmt::print(
      " Total Iterations                  = {}\n", settings::n_batches);
    fmt::print(
      " Number of Rays per Iteration      = {}\n", settings::n_particles);
    fmt::print(" Inactive Distance                 = {} cm\n",
      RandomRay::distance_inactive_);
    fmt::print(" Active Distance                   = {} cm\n",
      RandomRay::distance_active_);
    fmt::print(" Source Regions (SRs)              = {}\n", n_source_regions);
    fmt::print(
      " SRs Containing External Sources   = {}\n", n_external_source_regions);
    fmt::print(" Total Geometric Intersections     = {:.4e}\n",
      static_cast<double>(total_geometric_intersections));
    fmt::print("   Avg per Iteration               = {:.4e}\n",
      static_cast<double>(total_geometric_intersections) / settings::n_batches);
    fmt::print("   Avg per Iteration per SR        = {:.2f}\n",
      static_cast<double>(total_geometric_intersections) /
        static_cast<double>(settings::n_batches) / n_source_regions);
    fmt::print(" Avg SR Miss Rate per Iteration    = {:.4f}%\n", avg_miss_rate);
    fmt::print(" Energy Groups                     = {}\n", negroups);
    fmt::print(
      " Total Integrations                = {:.4e}\n", total_integrations);
    fmt::print("   Avg per Iteration               = {:.4e}\n",
      total_integrations / settings::n_batches);

    if (mpi::n_procs > 1){
      fmt::print(" MPI Ranks                         = {}\n", mpi::n_procs);
      fmt::print(" Avg Number of Subdomain Crossings = {}\n", avg_num_comms);
      fmt::print(" Load per MPI Rank                 = Rank {}: {:.2f}%\n", 0, rank_load_[0]);
      for (int i = 1; i < mpi::n_procs; i++) {
        fmt::print("                                     Rank {}: {:.2f}%\n", i, rank_load_[i]);
      }
    }

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

    header("Timing Statistics", 4);
    show_time("Total time for initialization", time_initialize.elapsed());
    show_time("Reading cross sections", time_read_xs.elapsed(), 1);
    show_time("Total simulation time", time_total.elapsed());
    show_time("Transport sweep only", time_transport.elapsed(), 1);
    show_time("Source update only", time_update_src.elapsed(), 1);
    show_time("Tally conversion only", time_tallies.elapsed(), 1);
    show_time("Other iteration routines", misc_time, 1);
    if (settings::run_mode == RunMode::EIGENVALUE) {
      show_time("Time in inactive batches", time_inactive.elapsed());
    }
    show_time("Time in active batches", time_active.elapsed());
    show_time("Time writing statepoints", time_statepoint.elapsed());
    show_time("Total time for finalization", time_finalize.elapsed());
    show_time("Time per integration", time_per_integration);
    if (mpi::n_procs > 1){
      show_time("Time in ray communication", time_ray_comms.elapsed());
      show_time("Time in decomposition handling", time_decomposition_handling.elapsed() 
      + time_ray_buffering.elapsed());
    }
  }

  if (settings::verbosity >= 4 && settings::run_mode == RunMode::EIGENVALUE) {
    header("Results", 4);
    fmt::print(" k-effective                       = {:.5f} +/- {:.5f}\n",
      simulation::keff, simulation::keff_std);
  }
}

void RandomRaySimulation::transport_sweep() {

  // Start timer for transport
  simulation::time_transport.start();

  // Transport sweep over all random rays for the iteration
  #pragma omp parallel for schedule(dynamic)                                     \
    reduction(+ : total_geometric_intersections_)
        for (int i = 0; i < settings::n_particles; i++) {
          RandomRay ray(i, domain_.get());
          total_geometric_intersections_ += ray.transport_history_based_single_ray();
        }

  simulation::time_transport.stop();

}

void RandomRaySimulation::transport_sweep_decomp(RayBank& RB) {

  // Transport sweep over all random rays for the iteration
  #pragma omp parallel for schedule(static)
    for (int i = 0; i < simulation::work_per_rank; i++) {
      uint64_t id = simulation::work_index[mpi::rank] + i;
      RandomRay ray(id, domain_.get());

      // Add ray to ray bank 
      #pragma omp critical (raybank)
      {
        RB.add_ray_to_bank(ray);
      }
    }

    // printf("RANK %d: Ray bank contains %d rays \n", mpi::rank, RB.ray_bank_size());

    int num_comms = 0;
    bool keep_running = true;

    // if (simulation::current_batch == 2){
    //   printf("Rank %d Starting transport sweep with %d rays \n", mpi::rank, RB.ray_bank_size());
    // }
    // printf("Rank %d Test2 \n", mpi::rank);

    // Move rays across ranks until they are terminated
    while (keep_running) {

      // Start timer for transport
      simulation::time_transport.start();

      // printf("RANK %d: Starting transport sweep with %d rays \n", mpi::rank, RB.ray_bank_size());

      #pragma omp parallel for schedule(dynamic)                                     \
        reduction(+ : total_geometric_intersections_)
          for (int i = 0; i < RB.ray_bank_size(); i++) {
            RandomRay& ray = RB.my_ray_list_[i];
            total_geometric_intersections_ += ray.transport_history_based_single_ray();
            // printf("RANK %d: Transported ray %ld \n", mpi::rank, ray.id());

            // printf("RANK %d: Transport of ray %d complete, new owner: %d \n", mpi::rank, i, ray.owner_rank_);


            // If ray has left my subdomain, buffer ray state
            if(ray.has_left_subdomain()){
              #pragma omp critical (raybuffer)
              {
                // printf("RANK %d: Buffering ray %d, new owner: %d \n", mpi::rank, i, ray.owner_rank_);
                // printf("Rank %d Test2.2 \n", mpi::rank);
                simulation::time_ray_buffering.start();
                num_ray_crossings_ += 1;
                RB.buffer_ray_data_to_send(ray, domain_.get());
                simulation::time_ray_buffering.stop();
                // printf("Rank %d Test2.3 \n", mpi::rank);
              }
            }
          }
      simulation::time_transport.stop();

      // Remove dead rays and move rays across subdomains between ranks
      simulation::time_decomposition_handling.start();

      // Check if new source regions were discovered and add them to subdomain map
      // if (mpi::decomp_map.new_sr_discovered()){
      //   mpi::decomp_map.exchange_sr_info()
      // }

      // printf("Rank %d Test3 \n", mpi::rank);

      // if (simulation::current_batch == 2){
      //   printf("Rank %d Completed transport of %d rays \n", mpi::rank, RB.ray_bank_size());
      // }

      // printf("RANK %d: Update ray bank \n", mpi::rank);
      RB.update(domain_.get());

      // bool test  = domain_->discovered_source_regions_.contains(SourceRegionKey(4, 104171));
      // printf("1 RANK %d: Source region (4, 104171) discovered? %d \n", mpi::rank, test);

      // bool test2 = false;
      // if(domain_->source_region_map_.find(SourceRegionKey(4, 104171)) != domain_->source_region_map_.end()){
      //   test2 = true;
      // }
      // printf("1 RANK %d: Source region (4, 104171) in main map? %d \n", mpi::rank, test2);


      // if (simulation::current_batch == 2){
      //   printf("Rank %d Completed send/recv, now has %d rays \n", mpi::rank, RB.ray_bank_size());
      // }

      // printf("Rank %d Test4 \n", mpi::rank);
      
      keep_running = RB.is_any_ray_alive();
      num_comms ++;
      simulation::time_decomposition_handling.stop();
    }

    // printf("Rank %d Test5 \n", mpi::rank);

    simulation::time_decomposition_handling.start();

    if (mpi::master) {
      printf("Max subdomain crossings: %d\n", num_comms);
      printf("Number of ray crossings: %ld\n", num_ray_crossings_);
    }

    // Calculate load per rank based on geometric intersections
    mpi::decomp_map.calculate_rank_load(total_geometric_intersections_);
    // mpi::decomp_map.calculate_rank_load(domain_.get());

    avg_num_comms_ += num_comms;

    // Reset ray bank list for next batch
    RB.reset_my_ray_list();

    simulation::time_decomposition_handling.stop();
}

} // namespace openmc
