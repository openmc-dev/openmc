#include "openmc/device_alloc.h"

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

} // namespace openmc
