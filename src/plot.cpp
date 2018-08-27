#include "openmc/plot.h"

namespace openmc {

void
voxel_init(hid_t file_id, const hsize_t* dims, hid_t* dspace, hid_t* dset,
           hid_t* memspace)
{
  // Create dataspace/dataset for voxel data
  *dspace = H5Screate_simple(3, dims, nullptr);
  *dset = H5Dcreate(file_id, "data", H5T_NATIVE_INT, *dspace, H5P_DEFAULT,
                    H5P_DEFAULT, H5P_DEFAULT);

  // Create dataspace for a slice of the voxel
  hsize_t dims_slice[2] {dims[1], dims[2]};
  *memspace = H5Screate_simple(2, dims_slice, nullptr);

  // Select hyperslab in dataspace
  hsize_t start[3] {0, 0, 0};
  hsize_t count[3] {1, dims[1], dims[2]};
  H5Sselect_hyperslab(*dspace, H5S_SELECT_SET, start, nullptr, count, nullptr);
}


void
voxel_write_slice(int x, hid_t dspace, hid_t dset, hid_t memspace, void* buf)
{
  hssize_t offset[3] {x - 1, 0, 0};
  H5Soffset_simple(dspace, offset);
  H5Dwrite(dset, H5T_NATIVE_INT, memspace, dspace, H5P_DEFAULT, buf);
}


void
voxel_finalize(hid_t dspace, hid_t dset, hid_t memspace)
{
  H5Dclose(dset);
  H5Sclose(dspace);
  H5Sclose(memspace);
}

} // namespace openmc
