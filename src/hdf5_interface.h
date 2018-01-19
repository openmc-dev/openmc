#ifndef HDF5_INTERFACE_H
#define HDF5_INTERFACE_H

#include <array>  // For std::array
#include <string.h>  // For strlen

#include "hdf5.h"


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


void
write_string(hid_t group_id, char const *name, char const *buffer)
{
  size_t buffer_len = strlen(buffer);
  hid_t datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, buffer_len);

  hid_t dataspace = H5Screate(H5S_SCALAR);

  hid_t dataset = H5Dcreate(group_id, name, datatype, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  H5Tclose(datatype);
  H5Sclose(dataspace);
  H5Dclose(dataset);
}

#endif //HDF5_INTERFACE_H
