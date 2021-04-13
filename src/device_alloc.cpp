#include "openmc/device_alloc.h"

#include "openmc/cell.h"
#include "openmc/lattice.h"
#include "openmc/nuclide.h"

namespace openmc {

  /*
void * device_alloc( size_t sz, int device_id )
{
  void * ptr = NULL;

  if( sz > 0 )
  {
    ptr = (void *) malloc(sz);
    assert(ptr != NULL);
  }

  return ptr;
}

void device_memcpy( void * dst_ptr, void * src_ptr, size_t sz, int dst_id, int src_id)
{
  memcpy(dst_ptr, src_ptr, sz);
}
*/

void move_read_only_data_to_device()
{
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
