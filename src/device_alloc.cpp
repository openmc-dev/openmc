#include "openmc/device_alloc.h"
#include "openmc/cell.h"
#include "openmc/lattice.h"

namespace openmc {

void * device_alloc( size_t sz, int device_id )
{
  void * ptr = NULL;

  if( sz > 0 )
  {
    #ifdef USE_DEVICE
    ptr = (void *) omp_target_alloc(sz, device_id);
    #else
    ptr = (void *) malloc(sz);
    #endif
    assert(ptr != NULL);
  }

  return ptr;
}

void device_memcpy( void * dst_ptr, void * src_ptr, size_t sz, int dst_id, int src_id)
{
  #ifdef USE_DEVICE
  omp_target_memcpy(dst_ptr, src_ptr, sz, 0, 0, dst_id, src_id);
  #else
  memcpy(dst_ptr, src_ptr, sz);
  #endif
}
  
void move_read_only_data_to_device(void)
{
  int host_id = omp_get_initial_device();
  int device_id = omp_get_default_device();
  size_t sz;
  
  // Surfaces ////////////////////////////////////////////////////////

  // Allocate space for all surfaces on device
  std::cout << "Moving " << model::surfaces.size() << " surfaces to device..." << std::endl;
  sz = model::surfaces.size() * sizeof(Surface);
  model::device_surfaces = (Surface *) device_alloc(sz, device_id);
  device_memcpy(model::device_surfaces, model::surfaces.data(), sz, device_id, host_id);
  
  // Universes ///////////////////////////////////////////////////////

  // Allocate and move vectors internal to each cell on the device
  std::cout << "Moving " << model::universes.size() << " universes to device..." << std::endl;
  for( auto& universe : model::universes ) {
    universe.allocate_and_copy_to_device();
  }
  // Move all universe data to device
  sz = model::universes.size() * sizeof(Universe);
  model::device_universes = (Universe *) device_alloc(sz, device_id);
  device_memcpy(model::device_universes, model::universes.data(), sz, device_id, host_id);
  
  // Cells ///////////////////////////////////////////////////////////

  // Allocate and move vectors internal to each cell on the device
  std::cout << "Moving " << model::cells.size() << " cells to device..." << std::endl;
  for( auto& cell : model::cells ) {
    cell.allocate_on_device();
    cell.copy_to_device();
  }
  // Move all cell data to device
  sz = model::cells.size() * sizeof(Cell);
  model::device_cells = (Cell *) device_alloc(sz, device_id);
  device_memcpy(model::device_cells, model::cells.data(), sz, device_id, host_id);
  
  // Lattices /////////////////////////////////////////////////////////

  // Allocate and move vectors internal to each cell on the device
  std::cout << "Moving " << model::lattices.size() << " lattices to device..." << std::endl;
  for( auto& lattice : model::lattices ) {
    lattice.allocate_and_copy_to_device();
  }
  // Move all lattice data to device
  sz = model::lattices.size() * sizeof(Lattice);
  model::device_lattices = (Lattice *) device_alloc(sz, device_id);
  device_memcpy(model::device_lattices, model::lattices.data(), sz, device_id, host_id);
}

} // namespace openmc
