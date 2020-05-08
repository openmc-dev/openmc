#ifndef OPENMC_DEVICE_ALLOC_H
#define OPENMC_DEVICE_ALLOC_H

#include <omp.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

//#define USE_DEVICE

namespace openmc {

void * device_alloc( size_t sz, int device_id );

void device_memcpy( void * dst_ptr, void * src_ptr, size_t sz, int dst_id, int src_id);

void move_read_only_data_to_device(void);

} // namespace openmc

#endif // OPENMC_DEVICE_ALLOC_H
