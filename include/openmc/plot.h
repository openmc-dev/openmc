#ifndef OPENMC_PLOT_H
#define OPENMC_PLOT_H

#include "hdf5.h"

namespace openmc {

extern "C" void voxel_init(hid_t file_id, const hsize_t* dims, hid_t* dspace,
                           hid_t* dset, hid_t* memspace);
extern "C" void voxel_write_slice(int x, hid_t dspace, hid_t dset,
                                  hid_t memspace, void* buf);
extern "C" void voxel_finalize(hid_t dspace, hid_t dset, hid_t memspace);

} // namespace openmc
#endif // OPENMC_PLOT_H
