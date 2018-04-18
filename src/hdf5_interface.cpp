#include "hdf5_interface.h"
#include "error.h"

#include "hdf5.h"
#include "hdf5_hl.h"

#include <array>
#include <string>
#include <sstream>

namespace openmc {

bool
using_mpio_device(hid_t obj_id)
{
  // Determine file that this object is part of
  hid_t file_id = H5Iget_file_id(obj_id);

  // Get file access property list
  hid_t fapl_id = H5Fget_access_plist(file_id);

  // Get low-level driver identifier
  hid_t driver = H5Pget_driver(fapl_id);

  // Free resources
  H5Pclose(fapl_id);
  H5Fclose(file_id);

  return driver == H5FD_MPIO;
}


hid_t
create_group(hid_t parent_id, char const *name)
{
  hid_t out = H5Gcreate(parent_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (out < 0) {
    std::stringstream err_msg;
    err_msg << "Failed to create HDF5 group \"" << name << "\"";
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
close_dataset(hid_t dataset_id)
{
  if (H5Dclose(dataset_id) < 0) fatal_error("Failed to close dataset");
}


void
close_group(hid_t group_id)
{
  if (H5Gclose(group_id) < 0) fatal_error("Failed to close group");
}


hid_t
file_open(const char* filename, char mode, bool parallel)
{
  bool create;
  unsigned int flags;
  switch (mode) {
  case 'r':
  case 'a':
    create = false;
    flags = (mode == 'r' ? H5F_ACC_RDONLY : H5F_ACC_RDWR);
    break;
  case 'w':
  case 'x':
    create = true;
    flags = (mode == 'x' ? H5F_ACC_EXCL : H5F_ACC_TRUNC);
    break;
  default:
    std::stringstream err_msg;
    err_msg <<  "Invalid file mode: " << mode;
    fatal_error(err_msg);
  }

  hid_t plist = H5P_DEFAULT;
#ifdef PHDF5
  if (parallel) {
    // Setup file access property list with parallel I/O access
    plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, mpi_intracomm, MPI_INFO_NULL);
  }
#endif

  // Open the file collectively
  hid_t file_id;
  if (create) {
    file_id = H5Fcreate(filename, flags, H5P_DEFAULT, plist);
  } else {
    file_id = H5Fopen(filename, flags, plist);
  }

#ifdef PHDF5
  // Close the property list
  if (parallel) H5Pclose(plist);
#endif

  return file_id;
}

hid_t
file_open(const std::string& filename, char mode, bool parallel=false)
{
  file_open(filename.c_str(), mode, parallel);
}

void file_close(hid_t file_id)
{
  H5Fclose(file_id);
}

bool
object_exists(hid_t object_id, const char* name)
{
  htri_t out = H5LTpath_valid(object_id, name, true);
  if (out < 0) {
    std::stringstream err_msg;
    err_msg << "Failed to check if object \"" << name << "\" exists.";
    fatal_error(err_msg);
  }
  return (out > 0);
}


hid_t
open_dataset(hid_t group_id, const char* name)
{
  if (object_exists(group_id, name)) {
    return H5Dopen(group_id, name, H5P_DEFAULT);
  } else {
    std::stringstream err_msg;
    err_msg << "Group \"" << name << "\" does not exist";
    fatal_error(err_msg);
  }
}


hid_t
open_group(hid_t group_id, const char* name)
{
  if (object_exists(group_id, name)) {
    return H5Gopen(group_id, name, H5P_DEFAULT);
  } else {
    std::stringstream err_msg;
    err_msg << "Group \"" << name << "\" does not exist";
    fatal_error(err_msg);
  }
}


void
read_double(hid_t obj_id, const char* name, double* buffer, bool indep)
{
  hid_t dset = obj_id;
  if (name) dset = H5Dopen(obj_id, name, H5P_DEFAULT);

  if (using_mpio_device(dset)) {
#ifdef PHDF5
    // Set up collective vs independent I/O
    auto data_xfer_mode {indep ? H5FD_MPIO_INDEPENDENT : H5FD_MPIO_COLLECTIVE};

    // Create dataset transfer property list
    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, data_xfer_mode);

    // Write data
    H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist, buffer);
    H5Pclose(plist);
#endif
  } else {
    H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  }

  if (name) H5Dclose(dset);
}


void
write_double(hid_t group_id, int ndim, const hsize_t* dims, const char* name,
             const double* buffer, bool indep)
{
  // If array is given, create a simple dataspace. Otherwise, create a scalar
  // datascape.
  hid_t dspace;
  if (ndim > 0) {
    dspace = H5Screate_simple(ndim, dims, nullptr);
  } else {
    dspace = H5Screate(H5S_SCALAR);
  }

  hid_t dset = H5Dcreate(group_id, name, H5T_NATIVE_DOUBLE, dspace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if (using_mpio_device(group_id)) {
#ifdef PHDF5
    // Set up collective vs independent I/O
    auto data_xfer_mode {indep ? H5FD_MPIO_INDEPENDENT : H5FD_MPIO_COLLECTIVE};

    // Create dataset transfer property list
    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, data_xfer_mode);

    // Write data
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist, buffer);
    H5Pclose(plist);
#endif
  } else {
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  }

  // Free resources
  H5Dclose(dset);
  H5Sclose(dspace);
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

}
