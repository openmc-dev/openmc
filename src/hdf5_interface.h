#ifndef HDF5_INTERFACE_H
#define HDF5_INTERFACE_H

#include <array>
#include <string.h>

#include "hdf5.h"

#include "error.h"


namespace openmc {


hid_t
create_group(hid_t parent_id, char const *name)
{
  hid_t out = H5Gcreate(parent_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (out < 0) {
    std::string err_msg{"Failed to create HDF5 group \""};
    err_msg += name;
    err_msg += "\"";
    fatal_error(err_msg);
  }
  return out;
}


hid_t
create_group(hid_t parent_id, const std::string &name)
{
  return create_group(parent_id, name.c_str());
}


void
close_group(hid_t group_id)
{
  herr_t err = H5Gclose(group_id);
  if (err < 0) {
    fatal_error("Failed to close HDF5 group");
  }
}


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


void
write_string(hid_t group_id, char const *name, const std::string &buffer)
{
  write_string(group_id, name, buffer.c_str());
}

} // namespace openmc
#endif //HDF5_INTERFACE_H
