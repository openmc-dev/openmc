#ifndef OPENMC_DEVICE_ALLOC_H
#define OPENMC_DEVICE_ALLOC_H

#include <omp.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

namespace openmc {

void move_read_only_data_to_device();

void release_data_from_device();

} // namespace openmc

#endif // OPENMC_DEVICE_ALLOC_H
