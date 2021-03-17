#include "openmc/device_alloc.h"
#include "openmc/cell.h"
#include "openmc/lattice.h"

namespace openmc {

void * device_alloc( size_t sz, int device_id )
{
  void * ptr = NULL;

  if( sz > 0 )
  {
    /*
    #ifdef USE_DEVICE
    ptr = (void *) omp_target_alloc(sz, device_id);
    #else
    ptr = (void *) malloc(sz);
    #endif
    */
    ptr = (void *) malloc(sz);
    assert(ptr != NULL);
  }

  return ptr;
}

void device_memcpy( void * dst_ptr, void * src_ptr, size_t sz, int dst_id, int src_id)
{
  /*
  #ifdef USE_DEVICE
  omp_target_memcpy(dst_ptr, src_ptr, sz, 0, 0, dst_id, src_id);
  #else
  memcpy(dst_ptr, src_ptr, sz);
  #endif
  */
  memcpy(dst_ptr, src_ptr, sz);
}
  
void move_read_only_data_to_device(void)
{
  int host_id = omp_get_initial_device();
  int device_id = omp_get_default_device();
  size_t sz;
  
  // Surfaces ////////////////////////////////////////////////////////

  // Allocate space for all surfaces on device
  //std::cout << "Moving " << model::surfaces.size() << " surfaces to device..." << std::endl;
  /*
  sz = model::surfaces.size() * sizeof(Surface);
  model::device_surfaces = (Surface *) device_alloc(sz, device_id);
  device_memcpy(model::device_surfaces, model::surfaces.data(), sz, device_id, host_id);
  */

  std::cout << "Moving " << model::surfaces.size() << " surfaces to device..." << std::endl;
  model::device_surfaces = model::surfaces.data();
  #pragma omp target enter data map(to: model::device_surfaces[:model::surfaces.size()])
  
  // Universes ///////////////////////////////////////////////////////
  
  model::device_universes = model::universes.data();
  #pragma omp target enter data map(to: model::device_universes[:model::universes.size()])

  // Allocate and move vectors internal to each cell on the device
  std::cout << "Moving " << model::universes.size() << " universes to device..." << std::endl;
  for( auto& universe : model::universes ) {
    universe.allocate_and_copy_to_device();
  }

  // Move all universe data to device
  /*
  sz = model::universes.size() * sizeof(Universe);
  model::device_universes = (Universe *) device_alloc(sz, device_id);
  device_memcpy(model::device_universes, model::universes.data(), sz, device_id, host_id);
  */
  
  // Cells ///////////////////////////////////////////////////////////
  /*
  Cell * tmp_cells = (Cell *) malloc(model::cells.size() * sizeof(Cell));
  memcpy(tmp_cells, model::cells.data(), model::cells.size() * sizeof(Cell));

  int * basic_ids = (int *) malloc( model::cells.size() * sizeof(int));
  for( int i =0; i < model::cells.size(); i++ )
    basic_ids[i] = tmp_cells[i].id_;

  //#pragma omp target enter data map(to: tmp_cells[:model::cells.size()])
    printf("Host   cell 0 ID: %d\n", tmp_cells[0].id_);
    #pragma omp target map(to: basic_ids[0:model::cells.size()]) map(to: tmp_cells[0:model::cells.size()])
    {
      printf("hello from the device! cell_id = %d\n", basic_ids[0]);
    }
  //#pragma omp target map(to: tmp_cells[:model::cells.size()])
  #pragma omp target map(to: tmp_cells[0:4])
  {
    printf("Device cell 0 ID: %d\n", tmp_cells[0].id_);
  }
  exit(0);
  */
  
  
  printf("moving %d global cells to device..\n", model::cells.size());
  model::device_cells = model::cells.data();
  #pragma omp target enter data map(to: model::device_cells[0:model::cells.size()])
  
  // Let's do a sanity check on our cells
  /*
  printf("Host    cell 0 ID: %d\n", model::device_cells[0].id_);
  #pragma omp target
  {
    //printf("Host    cell 0 ID: %d\n", model::device_cells[0].id_);
    model::device_cells[0].id_ *= 2;
  }
  */

  // Allocate and move vectors internal to each cell on the device
  //std::cout << "Moving " << model::cells.size() << " cells to device..." << std::endl;
  for( auto& cell : model::cells ) {
    //cell.allocate_on_device();
    cell.copy_to_device();
  }
  // Move all cell data to device
  /*
  sz = model::cells.size() * sizeof(Cell);
  model::device_cells = (Cell *) device_alloc(sz, device_id);
  device_memcpy(model::device_cells, model::cells.data(), sz, device_id, host_id);
  */
  
  // Lattices /////////////////////////////////////////////////////////

  model::device_lattices = model::lattices.data();
  #pragma omp target enter data map(to: model::device_lattices[:model::lattices.size()])

  // Allocate and move vectors internal to each cell on the device
  std::cout << "Moving " << model::lattices.size() << " lattices to device..." << std::endl;
  for( auto& lattice : model::lattices ) {
    lattice.allocate_and_copy_to_device();
  }
  // Move all lattice data to device
  /*
  sz = model::lattices.size() * sizeof(Lattice);
  model::device_lattices = (Lattice *) device_alloc(sz, device_id);
  device_memcpy(model::device_lattices, model::lattices.data(), sz, device_id, host_id);
  */
}

} // namespace openmc
