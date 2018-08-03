#include "hdf5_interface.h"

#include <array>
#include <cstring>
#include <sstream>
#include <string>

#include "hdf5.h"
#include "hdf5_hl.h"
#ifdef OPENMC_MPI
#include "mpi.h"
#include "message_passing.h"
#endif
#include "error.h"


namespace openmc {

bool
attribute_exists(hid_t obj_id, const char* name)
{
  htri_t out = H5Aexists_by_name(obj_id, ".", name, H5P_DEFAULT);
  return out > 0;
}


size_t
attribute_typesize(hid_t obj_id, const char* name)
{
  hid_t attr = H5Aopen(obj_id, name, H5P_DEFAULT);
  hid_t filetype = H5Aget_type(attr);
  size_t n = H5Tget_size(filetype);
  H5Tclose(filetype);
  H5Aclose(attr);
  return n;
}


void
get_shape(hid_t obj_id, hsize_t* dims)
{
  auto type = H5Iget_type(obj_id);
  hid_t dspace;
  if (type == H5I_DATASET) {
    dspace = H5Dget_space(obj_id);
  } else if (type == H5I_ATTR) {
    dspace = H5Aget_space(obj_id);
  }
  H5Sget_simple_extent_dims(dspace, dims, nullptr);
  H5Sclose(dspace);
}


void
get_shape_attr(hid_t obj_id, const char* name, hsize_t* dims)
{
  hid_t attr = H5Aopen(obj_id, name, H5P_DEFAULT);
  hid_t dspace = H5Aget_space(attr);
  H5Sget_simple_extent_dims(dspace, dims, nullptr);
  H5Sclose(dspace);
  H5Aclose(attr);
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


int
dataset_ndims(hid_t dset)
{
  hid_t dspace = H5Dget_space(dset);
  int ndims = H5Sget_simple_extent_ndims(dspace);
  H5Sclose(dspace);
  return ndims;
}


size_t
dataset_typesize(hid_t dset)
{
  hid_t filetype = H5Dget_type(dset);
  size_t n = H5Tget_size(filetype);
  H5Tclose(filetype);
  return n;
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
    H5Pset_fapl_mpio(plist, openmc::mpi::intracomm, MPI_INFO_NULL);
  }
#endif

  // Open the file collectively
  hid_t file_id;
  if (create) {
    file_id = H5Fcreate(filename, flags, H5P_DEFAULT, plist);
  } else {
    file_id = H5Fopen(filename, flags, plist);
  }
  if (file_id < 0) {
    std::stringstream msg;
    msg << "Failed to open HDF5 file with mode '" << mode << "': " << filename;
    fatal_error(msg);
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


void
get_name(hid_t obj_id, char* name)
{
  size_t size = 1 + H5Iget_name(obj_id, nullptr, 0);
  H5Iget_name(obj_id, name, size);
}


int get_num_datasets(hid_t group_id)
{
  // Determine number of links in the group
  H5G_info_t info;
  H5Gget_info(group_id, &info);

  // Iterate over links to get number of groups
  H5O_info_t oinfo;
  int ndatasets = 0;
  for (hsize_t i = 0; i < info.nlinks; ++i) {
    // Determine type of object (and skip non-group)
    H5Oget_info_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i, &oinfo,
                       H5P_DEFAULT);
    if (oinfo.type == H5O_TYPE_DATASET) ndatasets += 1;
  }

  return ndatasets;
}


int get_num_groups(hid_t group_id)
{
  // Determine number of links in the group
  H5G_info_t info;
  H5Gget_info(group_id, &info);

  // Iterate over links to get number of groups
  H5O_info_t oinfo;
  int ngroups = 0;
  for (hsize_t i = 0; i < info.nlinks; ++i) {
    // Determine type of object (and skip non-group)
    H5Oget_info_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i, &oinfo,
                       H5P_DEFAULT);
    if (oinfo.type == H5O_TYPE_GROUP) ngroups += 1;
  }

  return ngroups;
}


void
get_datasets(hid_t group_id, char* name[])
{
  // Determine number of links in the group
  H5G_info_t info;
  H5Gget_info(group_id, &info);

  // Iterate over links to get names
  H5O_info_t oinfo;
  hsize_t count = 0;
  size_t size;
  for (hsize_t i = 0; i < info.nlinks; ++i) {
    // Determine type of object (and skip non-group)
    H5Oget_info_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i, &oinfo,
                       H5P_DEFAULT);
    if (oinfo.type != H5O_TYPE_DATASET) continue;

    // Get size of name
    size = 1 + H5Lget_name_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC,
                                  i, nullptr, 0, H5P_DEFAULT);

    // Read name
    H5Lget_name_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i,
                       name[count], size, H5P_DEFAULT);
    count += 1;
  }
}


void
get_groups(hid_t group_id, char* name[])
{
  // Determine number of links in the group
  H5G_info_t info;
  H5Gget_info(group_id, &info);

  // Iterate over links to get names
  H5O_info_t oinfo;
  hsize_t count = 0;
  size_t size;
  for (hsize_t i = 0; i < info.nlinks; ++i) {
    // Determine type of object (and skip non-group)
    H5Oget_info_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i, &oinfo,
                       H5P_DEFAULT);
    if (oinfo.type != H5O_TYPE_GROUP) continue;

    // Get size of name
    size = 1 + H5Lget_name_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC,
                                  i, nullptr, 0, H5P_DEFAULT);

    // Read name
    H5Lget_name_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i,
                       name[count], size, H5P_DEFAULT);
    count += 1;
  }
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
read_attr(hid_t obj_id, const char* name, hid_t mem_type_id, void* buffer)
{
  hid_t attr = H5Aopen(obj_id, name, H5P_DEFAULT);
  H5Aread(attr, mem_type_id, buffer);
  H5Aclose(attr);
}


void
read_attr_double(hid_t obj_id, const char* name, double* buffer)
{
  read_attr(obj_id, name, H5T_NATIVE_DOUBLE, buffer);
}


void
read_attr_int(hid_t obj_id, const char* name, int* buffer)
{
  read_attr(obj_id, name, H5T_NATIVE_INT, buffer);
}


void
read_attr_string(hid_t obj_id, const char* name, size_t slen, char* buffer)
{
  // Create datatype for a string
  hid_t datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, slen + 1);

  // Read data into buffer
  read_attr(obj_id, name, datatype, buffer);

  // Free resources
  H5Tclose(datatype);
}


void
read_dataset(hid_t obj_id, const char* name, hid_t mem_type_id,
             void* buffer, bool indep)
{
  hid_t dset = obj_id;
  if (name) dset = open_dataset(obj_id, name);

  if (using_mpio_device(dset)) {
#ifdef PHDF5
    // Set up collective vs independent I/O
    auto data_xfer_mode = indep ? H5FD_MPIO_INDEPENDENT : H5FD_MPIO_COLLECTIVE;

    // Create dataset transfer property list
    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, data_xfer_mode);

    // Read data
    H5Dread(dset, mem_type_id, H5S_ALL, H5S_ALL, plist, buffer);
    H5Pclose(plist);
#endif
  } else {
    H5Dread(dset, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  }

  if (name) H5Dclose(dset);
}


void
read_double(hid_t obj_id, const char* name, double* buffer, bool indep)
{
  read_dataset(obj_id, name, H5T_NATIVE_DOUBLE, buffer, indep);
}


void
read_int(hid_t obj_id, const char* name, int* buffer, bool indep)
{
  read_dataset(obj_id, name, H5T_NATIVE_INT, buffer, indep);
}


void
read_llong(hid_t obj_id, const char* name, long long* buffer, bool indep)
{
  read_dataset(obj_id, name, H5T_NATIVE_LLONG, buffer, indep);
}


void
read_string(hid_t obj_id, const char* name, size_t slen, char* buffer, bool indep)
{
  // Create datatype for a string
  hid_t datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, slen + 1);

  // Read data into buffer
  read_dataset(obj_id, name, datatype, buffer, indep);

  // Free resources
  H5Tclose(datatype);
}


void
read_complex(hid_t obj_id, const char* name, double _Complex* buffer, bool indep)
{
  // Create compound datatype for complex numbers
  struct complex_t {
    double re;
    double im;
  };
  complex_t tmp;
  hid_t complex_id = H5Tcreate(H5T_COMPOUND, sizeof tmp);
  H5Tinsert(complex_id, "r", HOFFSET(complex_t, re), H5T_NATIVE_DOUBLE);
  H5Tinsert(complex_id, "i", HOFFSET(complex_t, im), H5T_NATIVE_DOUBLE);

  // Read data
  read_dataset(obj_id, name, complex_id, buffer, indep);

  // Free resources
  H5Tclose(complex_id);
}


void
read_nd_vector(hid_t obj_id, const char* name, std::vector<double>& result,
               bool must_have)
{
  if (object_exists(obj_id, name)) {
    read_double(obj_id, name, result.data(), true);
  } else if (must_have) {
    fatal_error(std::string("Must provide " + std::string(name) + "!"));
  }
}


void
read_nd_vector(hid_t obj_id, const char* name,
               std::vector<std::vector<double> >& result, bool must_have)
{
  if (object_exists(obj_id, name)) {
    int dim1 = result.size();
    int dim2 = result[0].size();
    double temp_arr[dim1 * dim2];
    read_double(obj_id, name, temp_arr, true);

    int temp_idx = 0;
    for (int i = 0; i < dim1; i++) {
      for (int j = 0; j < dim2; j++) {
        result[i][j] = temp_arr[temp_idx++];
      }
    }
  } else if (must_have) {
    fatal_error(std::string("Must provide " + std::string(name) + "!"));
  }
}


void
read_nd_vector(hid_t obj_id, const char* name,
               std::vector<std::vector<int> >& result, bool must_have)
{
  if (object_exists(obj_id, name)) {
    int dim1 = result.size();
    int dim2 = result[0].size();
    int temp_arr[dim1 * dim2];
    read_int(obj_id, name, temp_arr, true);

    int temp_idx = 0;
    for (int i = 0; i < dim1; i++) {
      for (int j = 0; j < dim2; j++) {
        result[i][j] = temp_arr[temp_idx++];
      }
    }
  } else if (must_have) {
    fatal_error(std::string("Must provide " + std::string(name) + "!"));
  }
}


void
read_nd_vector(hid_t obj_id, const char* name,
               std::vector<std::vector<std::vector<double> > >& result,
               bool must_have)
{
  if (object_exists(obj_id, name)) {
    int dim1 = result.size();
    int dim2 = result[0].size();
    int dim3 = result[0][0].size();
    double temp_arr[dim1 * dim2 * dim3];
    read_double(obj_id, name, temp_arr, true);

    int temp_idx = 0;
    for (int i = 0; i < dim1; i++) {
      for (int j = 0; j < dim2; j++) {
        for (int k = 0; k < dim3; k++) {
          result[i][j][k] = temp_arr[temp_idx++];
        }
      }
    }
  } else if (must_have) {
    fatal_error(std::string("Must provide " + std::string(name) + "!"));
  }
}

void
read_nd_vector(hid_t obj_id, const char* name,
               std::vector<std::vector<std::vector<int> > >& result,
               bool must_have)
{
  if (object_exists(obj_id, name)) {
    int dim1 = result.size();
    int dim2 = result[0].size();
    int dim3 = result[0][0].size();
    int temp_arr[dim1 * dim2 * dim3];
    read_int(obj_id, name, temp_arr, true);

    int temp_idx = 0;
    for (int i = 0; i < dim1; i++) {
      for (int j = 0; j < dim2; j++) {
        for (int k = 0; k < dim3; k++) {
          result[i][j][k] = temp_arr[temp_idx++];
        }
      }
    }
  } else if (must_have) {
    fatal_error(std::string("Must provide " + std::string(name) + "!"));
  }
}

void
read_nd_vector(hid_t obj_id, const char* name,
               std::vector<std::vector<std::vector<std::vector<double> > > >& result,
               bool must_have)
{
  if (object_exists(obj_id, name)) {
    int dim1 = result.size();
    int dim2 = result[0].size();
    int dim3 = result[0][0].size();
    int dim4 = result[0][0][0].size();
    double temp_arr[dim1 * dim2 * dim3 * dim4];
    read_double(obj_id, name, temp_arr, true);

    int temp_idx = 0;
    for (int i = 0; i < dim1; i++) {
      for (int j = 0; j < dim2; j++) {
        for (int k = 0; k < dim3; k++) {
          for (int l = 0; l < dim4; l++) {
            result[i][j][k][l] = temp_arr[temp_idx++];
          }
        }
      }
    }
  } else if (must_have) {
    fatal_error(std::string("Must provide " + std::string(name) + "!"));
  }
}

void
read_nd_vector(hid_t obj_id, const char* name,
               std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > >& result,
               bool must_have)
{
  if (object_exists(obj_id, name)) {
    int dim1 = result.size();
    int dim2 = result[0].size();
    int dim3 = result[0][0].size();
    int dim4 = result[0][0][0].size();
    int dim5 = result[0][0][0][0].size();
    double temp_arr[dim1 * dim2 * dim3 * dim4 * dim5];
    read_double(obj_id, name, temp_arr, true);

    int temp_idx = 0;
    for (int i = 0; i < dim1; i++) {
      for (int j = 0; j < dim2; j++) {
        for (int k = 0; k < dim3; k++) {
          for (int l = 0; l < dim4; l++) {
            for (int m = 0; m < dim5; m++) {
              result[i][j][k][l][m] = temp_arr[temp_idx++];
            }
          }
        }
      }
    }
  } else if (must_have) {
    fatal_error(std::string("Must provide " + std::string(name) + "!"));
  }
}


void
read_tally_results(hid_t group_id, hsize_t n_filter, hsize_t n_score, double* results)
{
  // Create dataspace for hyperslab in memory
  hsize_t dims[] {n_filter, n_score, 3};
  hsize_t start[] {0, 0, 1};
  hsize_t count[] {n_filter, n_score, 2};
  hid_t memspace = H5Screate_simple(3, dims, nullptr);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start, nullptr, count, nullptr);

  // Create and write dataset
  hid_t dset = H5Dopen(group_id, "results", H5P_DEFAULT);
  H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, H5S_ALL, H5P_DEFAULT, results);

  // Free resources
  H5Dclose(dset);
  H5Sclose(memspace);
}


void
write_attr(hid_t obj_id, int ndim, const hsize_t* dims, const char* name,
           hid_t mem_type_id, const void* buffer)
{
  // If array is given, create a simple dataspace. Otherwise, create a scalar
  // datascape.
  hid_t dspace;
  if (ndim > 0) {
    dspace = H5Screate_simple(ndim, dims, nullptr);
  } else {
    dspace = H5Screate(H5S_SCALAR);
  }

  // Create attribute and Write data
  hid_t attr = H5Acreate(obj_id, name, mem_type_id, dspace,
                         H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, mem_type_id, buffer);

  // Free resources
  H5Aclose(attr);
  H5Sclose(dspace);
}


void
write_attr_double(hid_t obj_id, int ndim, const hsize_t* dims, const char* name,
                  const double* buffer)
{
  write_attr(obj_id, ndim, dims, name, H5T_NATIVE_DOUBLE, buffer);
}


void
write_attr_int(hid_t obj_id, int ndim, const hsize_t* dims, const char* name,
               const int* buffer)
{
  write_attr(obj_id, ndim, dims, name, H5T_NATIVE_INT, buffer);
}


void
write_attr_string(hid_t obj_id, const char* name, const char* buffer)
{
  size_t n = strlen(buffer);
  if (n > 0) {
    // Set up appropriate datatype for a fixed-length string
    hid_t datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype, n);

    write_attr(obj_id, 0, nullptr, name, datatype, buffer);

    // Free resources
    H5Tclose(datatype);
  }
}


void
write_dataset(hid_t group_id, int ndim, const hsize_t* dims, const char* name,
              hid_t mem_type_id, const void* buffer, bool indep)
{
  // If array is given, create a simple dataspace. Otherwise, create a scalar
  // datascape.
  hid_t dspace;
  if (ndim > 0) {
    dspace = H5Screate_simple(ndim, dims, nullptr);
  } else {
    dspace = H5Screate(H5S_SCALAR);
  }

  hid_t dset = H5Dcreate(group_id, name, mem_type_id, dspace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if (using_mpio_device(group_id)) {
#ifdef PHDF5
    // Set up collective vs independent I/O
    auto data_xfer_mode = indep ? H5FD_MPIO_INDEPENDENT : H5FD_MPIO_COLLECTIVE;

    // Create dataset transfer property list
    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, data_xfer_mode);

    // Write data
    H5Dwrite(dset, mem_type_id, H5S_ALL, H5S_ALL, plist, buffer);
    H5Pclose(plist);
#endif
  } else {
    H5Dwrite(dset, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  }

  // Free resources
  H5Dclose(dset);
  H5Sclose(dspace);
}


void
write_double(hid_t group_id, int ndim, const hsize_t* dims, const char* name,
             const double* buffer, bool indep)
{
  write_dataset(group_id, ndim, dims, name, H5T_NATIVE_DOUBLE, buffer, indep);
}


void
write_int(hid_t group_id, int ndim, const hsize_t* dims, const char* name,
          const int* buffer, bool indep)
{
  write_dataset(group_id, ndim, dims, name, H5T_NATIVE_INT, buffer, indep);
}


void
write_llong(hid_t group_id, int ndim, const hsize_t* dims, const char* name,
            const long long* buffer, bool indep)
{
  write_dataset(group_id, ndim, dims, name, H5T_NATIVE_LLONG, buffer, indep);
}


void
write_string(hid_t group_id, int ndim, const hsize_t* dims, size_t slen,
             const char* name, const char* buffer, bool indep)
{
  if (slen > 0) {
    // Set up appropriate datatype for a fixed-length string
    hid_t datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype, slen);

    write_dataset(group_id, ndim, dims, name, datatype, buffer, indep);

    // Free resources
    H5Tclose(datatype);
  }
}


void
write_string(hid_t group_id, const char* name, const std::string& buffer, bool indep)
{
  write_string(group_id, 0, nullptr, buffer.length(), name, buffer.c_str(), indep);
}


void
write_tally_results(hid_t group_id, hsize_t n_filter, hsize_t n_score, const double* results)
{
  // Set dimensions of sum/sum_sq hyperslab to store
  hsize_t count[] {n_filter, n_score, 2};
  hid_t dspace = H5Screate_simple(3, count, nullptr);

  // Set dimensions of results array
  hsize_t dims[] {n_filter, n_score, 3};
  hsize_t start[] {0, 0, 1};
  hid_t memspace = H5Screate_simple(3, dims, nullptr);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start, nullptr, count, nullptr);

  // Create and write dataset
  hid_t dset = H5Dcreate(group_id, "results", H5T_NATIVE_DOUBLE, dspace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, H5S_ALL, H5P_DEFAULT, results);

  // Free resources
  H5Dclose(dset);
  H5Sclose(memspace);
  H5Sclose(dspace);
}


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

} // namespace openmc
