#include "openmc/device_alloc.h"

#include "openmc/bank.h"
#include "openmc/bremsstrahlung.h"
#include "openmc/cell.h"
#include "openmc/geometry.h"
#include "openmc/lattice.h"
#include "openmc/material.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/particle.h"
#include "openmc/photon.h"
#include "openmc/simulation.h"
#include "openmc/thermal.h"

#include "openmc/tallies/derivative.h"
#include "openmc/tallies/tally.h"


namespace openmc {

void enforce_assumptions()
{
  // TODO: These first two assumptions don't do anything, as no tallies are active at beginning of simulation

  // Notably, I have commented this capability out of particle::cross_vacuum_bc and particle::cross_reflective_bc
  assert(model::active_meshsurf_tallies.empty() && "Mesh surface tallies not yet supported.");

  // Commented out of particle::cross_reflective_bc
  assert(model::active_surface_tallies.empty() && "Surface tallies not yet supported.");

  // Assertions made when initializing particles
  assert(model::tally_derivs.size() <= FLUX_DERIVS_SIZE);
  assert(model::n_tally_filters <= FILTER_MATCHES_SIZE);
  assert(model::n_coord_levels <= COORD_SIZE);
  #ifndef NO_MICRO_XS_CACHE
  assert(data::nuclides_size <= NEUTRON_XS_SIZE);
  #endif
  assert(data::elements_size <= PHOTON_XS_SIZE);
}

void move_settings_to_device()
{
  // settings.h
  #pragma omp target update to(settings::dagmc)
  #pragma omp target update to(settings::run_CE)
  #pragma omp target update to(settings::max_lost_particles)
  #pragma omp target update to(settings::rel_max_lost_particles)
  #pragma omp target update to(settings::gen_per_batch)
  #pragma omp target update to(settings::run_mode)
  #pragma omp target update to(settings::n_particles)
  #pragma omp target update to(settings::temperature_method)
  #pragma omp target update to(settings::urr_ptables_on)
  #pragma omp target update to(settings::create_fission_neutrons)
  #pragma omp target update to(settings::survival_biasing)
  #pragma omp target update to(settings::res_scat_method)
  #pragma omp target update to(settings::res_scat_energy_min)
  #pragma omp target update to(settings::res_scat_energy_max)
  #pragma omp target update to(settings::weight_cutoff)
  #pragma omp target update to(settings::weight_survive)
  #pragma omp target update to(settings::electron_treatment)
  settings::energy_cutoff[0]; // Lazy extern template expansion workaround
  #pragma omp target update to(settings::energy_cutoff)
  #pragma omp target update to(settings::n_log_bins)
  #pragma omp target update to(settings::assume_separate)
  #pragma omp target update to(settings::check_overlaps)
  #pragma omp target update to(settings::max_particles_in_flight)
  #pragma omp target update to(settings::minimum_sort_items)

  // message_passing.h
  #pragma omp target update to(mpi::rank)
  #pragma omp target update to(mpi::n_procs)
  #pragma omp target update to(mpi::master)

  // simulation.h
  #pragma omp target update to(simulation::total_gen)
  #pragma omp target update to(simulation::current_batch)
  #pragma omp target update to(simulation::current_gen)
  #pragma omp target update to(simulation::total_weight)
  #pragma omp target update to(simulation::need_depletion_rx)
  #pragma omp target update to(simulation::log_spacing)

  // geometry.h
  #pragma omp target update to(model::root_universe)
}

void move_read_only_data_to_device()
{
  // Enforce any device-specific assumptions or limitations on user inputs
  enforce_assumptions();

  // Copy all global settings into device globals
  move_settings_to_device();

  int host_id = omp_get_initial_device();
  int device_id = omp_get_default_device();
  size_t sz;

  // Surfaces ////////////////////////////////////////////////////////

  std::cout << "Moving " << model::surfaces.size() << " surfaces to device..." << std::endl;
  model::device_surfaces = model::surfaces.data();
  #pragma omp target enter data map(to: model::device_surfaces[:model::surfaces.size()])

  // Universes ///////////////////////////////////////////////////////

  std::cout << "Moving " << model::universes.size() << " universes to device..." << std::endl;
  model::device_universes = model::universes.data();
  #pragma omp target enter data map(to: model::device_universes[:model::universes.size()])
  for( auto& universe : model::universes ) {
    universe.allocate_and_copy_to_device();
  }

  // Cells //////////////////////////////////////////////////////////

  std::cout << "Moving " << model::cells.size() << " cells to device..." << std::endl;
  model::device_cells = model::cells.data();
  #pragma omp target enter data map(to: model::device_cells[0:model::cells.size()])
  for( auto& cell : model::cells ) {
    cell.copy_to_device();
  }

  // Lattices /////////////////////////////////////////////////////////

  std::cout << "Moving " << model::lattices.size() << " lattices to device..." << std::endl;
  model::device_lattices = model::lattices.data();
  #pragma omp target enter data map(to: model::device_lattices[:model::lattices.size()])
  for( auto& lattice : model::lattices ) {
    lattice.allocate_and_copy_to_device();
  }

  // Nuclear data /////////////////////////////////////////////////////
  data::energy_min[0]; // Lazy extern template expansion workaround
  data::energy_max[0]; // Lazy extern template expansion workaround
  #pragma omp target update to(data::energy_min)
  #pragma omp target update to(data::energy_max)
  #pragma omp target update to(data::nuclides_size)

  // Flatten nuclides before copying
  for (int i = 0; i < data::nuclides_size; ++i) {
    auto& nuc = data::nuclides[i];

    // URR data flattening
    for (auto& u : nuc.urr_data_) {
      u.flatten_urr_data();
    }

    // Pointwise XS data flattening
    nuc.flatten_xs_data();

    // Windowed multipole
    nuc.flatten_wmp_data();
  }

  std::cout << "Moving " << data::nuclides_size << " nuclides to device..." << std::endl;
  #pragma omp target enter data map(to: data::nuclides[:data::nuclides_size])
  for (int i = 0; i < data::nuclides_size; ++i) {
    auto& nuc = data::nuclides[i];
    nuc.copy_to_device();
  }

  data::device_thermal_scatt = data::thermal_scatt.data();
  #pragma omp target enter data map(to: data::device_thermal_scatt[:data::thermal_scatt.size()])
  for (auto& ts : data::thermal_scatt) {
    ts.copy_to_device();
  }

  // Photon data /////////////////////////////////////////////////////

  #pragma omp target update to(data::compton_profile_pz_size)
  #pragma omp target enter data map(to: data::compton_profile_pz[:data::compton_profile_pz_size])

  #pragma omp target update to(data::elements_size)
  #pragma omp target enter data map(to: data::elements[:data::elements_size])
  for (int i = 0; i < data::elements_size; ++i) {
    auto& elm = data::elements[i];
    std::cout << "Moving " << elm.name_ << " data to device..." << std::endl;
    elm.copy_to_device();
  }
  data::device_ttb_e_grid = data::ttb_e_grid.data();
  #pragma omp target update to(data::ttb_e_grid_size)
  #pragma omp target enter data map(to: data::device_ttb_e_grid[:data::ttb_e_grid.size()])

  // Materials /////////////////////////////////////////////////////////

  std::cout << "Moving " << model::materials_size << " materials to device..." << std::endl;
  int min = 99999;
  int max = 0;
  int n_over_200 = 0;
  int n_under_200 = 0;
  #pragma omp target update to(model::materials_size)
  #pragma omp target enter data map(to: model::materials[:model::materials_size])
  for (int i = 0; i < model::materials_size; i++) {
    model::materials[i].copy_to_device();
    if(model::materials[i].fissionable())
    {
      int num_nucs = model::materials[i].nuclide_.size();
      if( num_nucs < min )
        min = num_nucs;
      if( num_nucs > max )
       max = num_nucs;
     if( num_nucs >= 200 )
       n_over_200++;
     else
       n_under_200++;
    }
  }
  std::cout << "Fissionable Material Statistics:" << std::endl <<
  " Max Nuclide Count: " << max << std::endl <<
  " Min Nuclide Count: " << min << std::endl <<
  " # Fissionable Materials with >= 200 Nuclides: " << n_over_200 << std::endl <<
  " # Fissionable Materials with  < 200 Nuclides: " << n_under_200 << std::endl;

  // Source Bank ///////////////////////////////////////////////////////

  simulation::device_source_bank = simulation::source_bank.data();
  #pragma omp target enter data map(alloc: simulation::device_source_bank[:simulation::source_bank.size()])
  simulation::fission_bank.allocate_on_device();

  // MPI Work Indices ///////////////////////////////////////////////////

  simulation::device_work_index = simulation::work_index.data();
  #pragma omp target enter data map(to: simulation::device_work_index[:simulation::work_index.size()])

  // Progeny per Particle ///////////////////////////////////////////////////

  simulation::device_progeny_per_particle = simulation::progeny_per_particle.data();
  #pragma omp target enter data map(alloc: simulation::device_progeny_per_particle[:simulation::progeny_per_particle.size()])

  // Filters ////////////////////////////////////////////////////////////////

  std::cout << "Moving " << model::n_tally_filters << " tally filters to device..." << std::endl;
  #pragma omp target update to(model::n_tally_filters)
  #pragma omp target enter data map(to: model::tally_filters[:model::n_tally_filters])
  for (int i = 0; i < model::n_tally_filters; i++) {
    model::tally_filters[i].copy_to_device();
  }

  // Meshes ////////////////////////////////////////////////////////////////

  std::cout << "Moving " << model::meshes_size << " meshes to device..." << std::endl;
  #pragma omp target update to(model::meshes_size)
  #pragma omp target enter data map(to: model::meshes[:model::meshes_size])
  for (int i = 0; i < model::meshes_size; i++) {
    model::meshes[i].copy_to_device();
  }

  // Tallies ///////////////////////////////////////////////////

  std::cout << "Moving " << model::tallies_size << " tallies to device..." << std::endl;
  #pragma omp target update to(model::tallies_size)
  #pragma omp target enter data map(to: model::tallies[:model::tallies_size])
  for (int i = 0; i < model::tallies_size; ++i) {
    auto& tally = model::tallies[i];
    std::cout << "Moving " << tally.name_ << " data to device..." << std::endl;
    tally.copy_to_device();
  }
  model::device_active_tallies = model::active_tallies.data();
  model::device_active_collision_tallies = model::active_collision_tallies.data();
  model::device_active_tracklength_tallies = model::active_tracklength_tallies.data();
  #pragma omp target enter data map(to: model::device_active_tallies[:model::active_tallies.size()])
  #pragma omp target enter data map(to: model::device_active_collision_tallies[:model::active_collision_tallies.size()])
  #pragma omp target enter data map(to: model::device_active_tracklength_tallies[:model::active_tracklength_tallies.size()])

  #ifdef OPENMC_MPI
  MPI_Barrier( mpi::intracomm );
  #endif
}


void release_data_from_device()
{
  std::cout << "Releasing data from device..." << std::endl;
  for (int i = 0; i < data::nuclides_size; ++i) {
    data::nuclides[i].release_from_device();
  }

  for (int i = 0; i < data::elements_size; ++i) {
    data::elements[i].release_from_device();
  }

  for (int i = 0; i < model::tallies_size; ++i) {
    model::tallies[i].release_from_device();
  }
}


} // namespace openmc
