#include "openmc/device_alloc.h"

#include "openmc/cell.h"
#include "openmc/lattice.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"

#include "openmc/tallies/tally.h"

namespace openmc {

void enforce_assumptions()
{
  // Notably, I have commented this capability out of particle::cross_vacuum_bc and particle::cross_reflective_bc
  assert(model::active_meshsurf_tallies.empty() && "Mesh surface tallies not yet supported.");
  
  // Commented out of particle::cross_reflective_bc
  assert(model::active_surface_tallies.empty() && "Surface tallies not yet supported.");
}

void move_settings_to_device()
{
  #pragma omp target update to(settings::dagmc)
  #pragma omp target update to(settings::run_CE)
  #pragma omp target update to(settings::max_lost_particles)
  #pragma omp target update to(settings::rel_max_lost_particles)
  #pragma omp target update to(settings::gen_per_batch)
  #pragma omp target update to(settings::run_mode)
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

  // Nuclides /////////////////////////////////////////////////////////
  for (auto& nuc : data::nuclides) {
    std::cout << "Moving " << nuc->name_ << " data to device..." << std::endl;
    for (auto& rx : nuc->reactions_) {
      for (auto& product : rx->products_) {
        for (auto& d : product.distribution_) {
          #pragma omp target enter data map(to: d)
          d.copy_to_device();
        }
      }
    }
  }
  
  // Materials /////////////////////////////////////////////////////////
  
  std::cout << "Moving " << model::materials.size() << " materials to device..." << std::endl;
  model::device_materials = model::materials.data();
  #pragma omp target enter data map(to: model::device_materials[:model::materials.size()])
  // TODO: Deep copy of materials
}

void release_data_from_device()
{
  std::cout << "Releasing data from device..." << std::endl;
  for (auto& nuc : data::nuclides) {
    for (auto& rx : nuc->reactions_) {
      for (auto& product : rx->products_) {
        for (auto& d : product.distribution_) {
          d.release_device();
          #pragma omp target exit data map(release: d)
        }
      }
    }
  }
}


} // namespace openmc
