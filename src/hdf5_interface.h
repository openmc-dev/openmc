#ifndef HDF5_INTERFACE_H
#define HDF5_INTERFACE_H

#include "hdf5.h"
#include "hdf5_hl.h"

#include <array>
#include <string>
#include <sstream>


namespace openmc {

bool using_mpio_device(hid_t obj_id);
hid_t create_group(hid_t parent_id, const char* name);
hid_t create_group(hid_t parent_id, const std::string& name);
void close_dataset(hid_t dataset_id);
void close_group(hid_t group_id);
extern "C" hid_t file_open(const char* filename, char mode, bool parallel);
hid_t file_open(const std::string& filename, char mode, bool parallel);
void file_close(hid_t file_id);
bool object_exists(hid_t object_id, const char* name);
hid_t open_dataset(hid_t group_id, const char* name);
hid_t open_group(hid_t group_id, const char* name);


template<std::size_t array_len> void
write_double_1D(hid_t group_id, char const *name,
                std::array<double, array_len> &buffer)
{
  hsize_t dims[1]{array_len};
  hid_t dataspace = H5Screate_simple(1, dims, NULL);

  hid_t dataset = H5Dcreate(group_id, name, H5T_NATIVE_DOUBLE, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           &buffer[0]);

  H5Sclose(dataspace);
  H5Dclose(dataset);
}

void write_double(hid_t group_id, int ndim, const hsize_t* dims, const char* name,
                  double* buffer, bool indep);

void write_string(hid_t group_id, char const *name, char const *buffer);
void write_string(hid_t group_id, char const *name, const std::string &buffer);

} // namespace openmc
#endif //HDF5_INTERFACE_H
