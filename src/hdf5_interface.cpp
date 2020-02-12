#include "openmc/hdf5_interface.h"

#include <array>
#include <cstring>
#include <stdexcept>
#include <string>

#include "xtensor/xtensor.hpp"
#include "xtensor/xarray.hpp"
#include <fmt/core.h>

#include "hdf5.h"
#include "hdf5_hl.h"
#ifdef OPENMC_MPI
#include "mpi.h"
#include "openmc/message_passing.h"
#endif


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
  } else {
    throw std::runtime_error{
      "Expected dataset or attribute in call to get_shape."};
  }
  H5Sget_simple_extent_dims(dspace, dims, nullptr);
  H5Sclose(dspace);
}


std::vector<hsize_t> attribute_shape(hid_t obj_id, const char* name)
{
  hid_t attr = H5Aopen(obj_id, name, H5P_DEFAULT);
  std::vector<hsize_t> shape = object_shape(attr);
  H5Aclose(attr);
  return shape;
}

std::vector<hsize_t> object_shape(hid_t obj_id)
{
  // Get number of dimensions
  auto type = H5Iget_type(obj_id);
  hid_t dspace;
  if (type == H5I_DATASET) {
    dspace = H5Dget_space(obj_id);
  } else if (type == H5I_ATTR) {
    dspace = H5Aget_space(obj_id);
  } else {
    throw std::runtime_error{
      "Expected dataset or attribute in call to object_shape."};
  }
  int n = H5Sget_simple_extent_ndims(dspace);

  // Get shape of array
  std::vector<hsize_t> shape(n);
  H5Sget_simple_extent_dims(dspace, shape.data(), nullptr);

  // Free resources and return
  H5Sclose(dspace);
  return shape;
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
    fatal_error(fmt::format("Failed to create HDF5 group \"{}\"", name));
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
dataset_typesize(hid_t obj_id, const char* name)
{
  hid_t dset = open_dataset(obj_id, name);
  hid_t filetype = H5Dget_type(dset);
  size_t n = H5Tget_size(filetype);
  H5Tclose(filetype);
  close_dataset(dset);
  return n;
}


void
ensure_exists(hid_t obj_id, const char* name, bool attribute)
{
  if (attribute) {
    if (!attribute_exists(obj_id, name)) {
      fatal_error(fmt::format("Attribute \"{}\" does not exist in object {}",
        name, object_name(obj_id)));
    }
  } else {
    if (!object_exists(obj_id, name)) {
      fatal_error(fmt::format("Object \"{}\" does not exist in object {}",
        name, object_name(obj_id)));
    }
  }
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
      fatal_error(fmt::format("Invalid file mode: ", mode));
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
    fatal_error(fmt::format(
      "Failed to open HDF5 file with mode '{}': {}", mode, filename));
  }

#ifdef PHDF5
  // Close the property list
  if (parallel) H5Pclose(plist);
#endif

  return file_id;
}

hid_t
file_open(const std::string& filename, char mode, bool parallel)
{
  return file_open(filename.c_str(), mode, parallel);
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

std::vector<std::string>
member_names(hid_t group_id, H5O_type_t type)
{
  // Determine number of links in the group
  H5G_info_t info;
  H5Gget_info(group_id, &info);

  // Iterate over links to get names
  H5O_info_t oinfo;
  size_t size;
  std::vector<std::string> names;
  for (hsize_t i = 0; i < info.nlinks; ++i) {
    // Determine type of object (and skip non-group)
    H5Oget_info_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i, &oinfo,
                       H5P_DEFAULT);
    if (oinfo.type != type) continue;

    // Get size of name
    size = 1 + H5Lget_name_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC,
                                  i, nullptr, 0, H5P_DEFAULT);

    // Read name
    char* buffer = new char[size];
    H5Lget_name_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i,
                       buffer, size, H5P_DEFAULT);
    names.emplace_back(&buffer[0]);
    delete[] buffer;
  }
  return names;
}

std::vector<std::string>
group_names(hid_t group_id)
{
  return member_names(group_id, H5O_TYPE_GROUP);
}

std::vector<std::string>
dataset_names(hid_t group_id)
{
  return member_names(group_id, H5O_TYPE_DATASET);
}

bool
object_exists(hid_t object_id, const char* name)
{
  htri_t out = H5LTpath_valid(object_id, name, true);
  if (out < 0) {
    fatal_error(fmt::format("Failed to check if object \"{}\" exists.", name));
  }
  return (out > 0);
}


std::string
object_name(hid_t obj_id)
{
  // Determine size and create buffer
  size_t size = 1 + H5Iget_name(obj_id, nullptr, 0);
  char* buffer = new char[size];

  // Read and return name
  H5Iget_name(obj_id, buffer, size);
  std::string str = buffer;
  delete[] buffer;
  return str;
}


hid_t
open_dataset(hid_t group_id, const char* name)
{
  ensure_exists(group_id, name);
  return H5Dopen(group_id, name, H5P_DEFAULT);
}


hid_t
open_group(hid_t group_id, const char* name)
{
  ensure_exists(group_id, name);
  return H5Gopen(group_id, name, H5P_DEFAULT);
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
  H5Tset_size(datatype, slen);
  // numpy uses null-padding when writing fixed-length strings
  H5Tset_strpad(datatype, H5T_STR_NULLPAD);

  // Read data into buffer
  read_attr(obj_id, name, datatype, buffer);

  // Free resources
  H5Tclose(datatype);
}


void
read_dataset_lowlevel(hid_t obj_id, const char* name, hid_t mem_type_id,
                      hid_t mem_space_id, bool indep, void* buffer)
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
    H5Dread(dset, mem_type_id, mem_space_id, H5S_ALL, plist, buffer);
    H5Pclose(plist);
#endif
  } else {
    H5Dread(dset, mem_type_id, mem_space_id, H5S_ALL, H5P_DEFAULT, buffer);
  }

  if (name) H5Dclose(dset);
}

template<>
void read_dataset(hid_t dset, xt::xarray<std::complex<double>>& arr, bool indep)
{
  // Get shape of dataset
  std::vector<hsize_t> shape = object_shape(dset);

  // Allocate new array to read data into
  std::size_t size = 1;
  for (const auto x : shape)
    size *= x;
  std::vector<std::complex<double>> buffer(size);

  // Read data from attribute
  read_complex(dset, nullptr, buffer.data(), indep);

  // Adapt into xarray
  arr = xt::adapt(buffer, shape);
}


void
read_double(hid_t obj_id, const char* name, double* buffer, bool indep)
{
  read_dataset_lowlevel(obj_id, name, H5T_NATIVE_DOUBLE, H5S_ALL, indep,
                        buffer);
}


void
read_int(hid_t obj_id, const char* name, int* buffer, bool indep)
{
  read_dataset_lowlevel(obj_id, name, H5T_NATIVE_INT, H5S_ALL, indep, buffer);
}


void
read_llong(hid_t obj_id, const char* name, long long* buffer, bool indep)
{
  read_dataset_lowlevel(obj_id, name, H5T_NATIVE_LLONG, H5S_ALL, indep, buffer);
}


void
read_string(hid_t obj_id, const char* name, size_t slen, char* buffer,
            bool indep)
{
  // Create datatype for a string
  hid_t datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, slen);
  // numpy uses null-padding when writing fixed-length strings
  H5Tset_strpad(datatype, H5T_STR_NULLPAD);

  // Read data into buffer
  read_dataset_lowlevel(obj_id, name, datatype, H5S_ALL, indep, buffer);

  // Free resources
  H5Tclose(datatype);
}


void
read_complex(hid_t obj_id, const char* name, std::complex<double>* buffer,
             bool indep)
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
  read_dataset_lowlevel(obj_id, name, complex_id, H5S_ALL, indep, buffer);

  // Free resources
  H5Tclose(complex_id);
}


void
read_tally_results(hid_t group_id, hsize_t n_filter, hsize_t n_score,
                   double* results)
{
  // Create dataspace for hyperslab in memory
  constexpr int ndim = 3;
  hsize_t dims[ndim] {n_filter, n_score, 3};
  hsize_t start[ndim] {0, 0, 1};
  hsize_t count[ndim] {n_filter, n_score, 2};
  hid_t memspace = H5Screate_simple(ndim, dims, nullptr);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start, nullptr, count, nullptr);

  // Read the dataset
  read_dataset_lowlevel(group_id, "results", H5T_NATIVE_DOUBLE, memspace,
                        false, results);

  // Free resources
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
write_dataset_lowlevel(hid_t group_id, int ndim, const hsize_t* dims,
  const char* name, hid_t mem_type_id, hid_t mem_space_id, bool indep,
  const void* buffer)
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
    H5Dwrite(dset, mem_type_id, mem_space_id, H5S_ALL, plist, buffer);
    H5Pclose(plist);
#endif
  } else {
    H5Dwrite(dset, mem_type_id, mem_space_id, H5S_ALL, H5P_DEFAULT, buffer);
  }

  // Free resources
  H5Dclose(dset);
  H5Sclose(dspace);
}


void
write_double(hid_t group_id, int ndim, const hsize_t* dims, const char* name,
             const double* buffer, bool indep)
{
  write_dataset_lowlevel(group_id, ndim, dims, name, H5T_NATIVE_DOUBLE, H5S_ALL,
                         indep, buffer);
}


void
write_int(hid_t group_id, int ndim, const hsize_t* dims, const char* name,
          const int* buffer, bool indep)
{
  write_dataset_lowlevel(group_id, ndim, dims, name, H5T_NATIVE_INT, H5S_ALL,
                         indep, buffer);
}


void
write_llong(hid_t group_id, int ndim, const hsize_t* dims, const char* name,
            const long long* buffer, bool indep)
{
  write_dataset_lowlevel(group_id, ndim, dims, name, H5T_NATIVE_LLONG, H5S_ALL,
                         indep, buffer);
}


void
write_string(hid_t group_id, int ndim, const hsize_t* dims, size_t slen,
             const char* name, const char* buffer, bool indep)
{
  if (slen > 0) {
    // Set up appropriate datatype for a fixed-length string
    hid_t datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype, slen);

    write_dataset_lowlevel(group_id, ndim, dims, name, datatype, H5S_ALL, indep,
                           buffer);

    // Free resources
    H5Tclose(datatype);
  }
}


void
write_string(hid_t group_id, const char* name, const std::string& buffer,
             bool indep)
{
  write_string(group_id, 0, nullptr, buffer.length(), name, buffer.c_str(),
               indep);
}


void
write_tally_results(hid_t group_id, hsize_t n_filter, hsize_t n_score,
                    const double* results)
{
  // Set dimensions of sum/sum_sq hyperslab to store
  constexpr int ndim = 3;
  hsize_t count[ndim] {n_filter, n_score, 2};

  // Set dimensions of results array
  hsize_t dims[ndim] {n_filter, n_score, 3};
  hsize_t start[ndim] {0, 0, 1};
  hid_t memspace = H5Screate_simple(ndim, dims, nullptr);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start, nullptr, count, nullptr);

  // Create and write dataset
  write_dataset_lowlevel(group_id, ndim, count, "results", H5T_NATIVE_DOUBLE,
                         memspace, false, results);

  // Free resources
  H5Sclose(memspace);
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

// Specializations of the H5TypeMap template struct
template<>
const hid_t H5TypeMap<bool>::type_id = H5T_NATIVE_INT8;
template<>
const hid_t H5TypeMap<int>::type_id = H5T_NATIVE_INT;
template<>
const hid_t H5TypeMap<unsigned long>::type_id = H5T_NATIVE_ULONG;
template<>
const hid_t H5TypeMap<unsigned long long>::type_id = H5T_NATIVE_ULLONG;
template<>
const hid_t H5TypeMap<unsigned int>::type_id = H5T_NATIVE_UINT;
template<>
const hid_t H5TypeMap<int64_t>::type_id = H5T_NATIVE_INT64;
template<>
const hid_t H5TypeMap<double>::type_id = H5T_NATIVE_DOUBLE;
template <>
const hid_t H5TypeMap<char>::type_id = H5T_NATIVE_CHAR;

} // namespace openmc
