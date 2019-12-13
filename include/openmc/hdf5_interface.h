#ifndef OPENMC_HDF5_INTERFACE_H
#define OPENMC_HDF5_INTERFACE_H

#include <algorithm> // for min
#include <array>
#include <complex>
#include <cstddef>
#include <cstring> // for strlen
#include <string>
#include <sstream>
#include <type_traits>
#include <vector>

#include "hdf5.h"
#include "hdf5_hl.h"
#include "xtensor/xadapt.hpp"
#include "xtensor/xarray.hpp"

#include "openmc/position.h"
#include "openmc/error.h"

namespace openmc {

//==============================================================================
// Low-level internal functions
//==============================================================================

void read_attr(hid_t obj_id, const char* name, hid_t mem_type_id,
  void* buffer);

void write_attr(hid_t obj_id, int ndim, const hsize_t* dims, const char* name,
  hid_t mem_type_id, const void* buffer);

void read_dataset_lowlevel(hid_t obj_id, const char* name, hid_t mem_type_id,
  hid_t mem_space_id, bool indep, void* buffer);

void write_dataset_lowlevel(hid_t group_id, int ndim, const hsize_t* dims,
  const char* name, hid_t mem_type_id, hid_t mem_space_id, bool indep,
  const void* buffer);

bool using_mpio_device(hid_t obj_id);

//==============================================================================
// Normal functions that are used to read/write files
//==============================================================================

hid_t create_group(hid_t parent_id, const std::string& name);

inline hid_t create_group(hid_t parent_id, const std::stringstream& name)
{return create_group(parent_id, name.str());}


hid_t file_open(const std::string& filename, char mode, bool parallel=false);

void write_string(hid_t group_id, const char* name, const std::string& buffer,
                  bool indep);

std::vector<hsize_t> attribute_shape(hid_t obj_id, const char* name);
std::vector<std::string> dataset_names(hid_t group_id);
void ensure_exists(hid_t obj_id, const char* name, bool attribute=false);
std::vector<std::string> group_names(hid_t group_id);
std::vector<hsize_t> object_shape(hid_t obj_id);
std::string object_name(hid_t obj_id);

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  bool attribute_exists(hid_t obj_id, const char* name);
  size_t attribute_typesize(hid_t obj_id, const char* name);
  hid_t create_group(hid_t parent_id, const char* name);
  void close_dataset(hid_t dataset_id);
  void close_group(hid_t group_id);
  int dataset_ndims(hid_t dset);
  size_t dataset_typesize(hid_t obj_id, const char* name);
  hid_t file_open(const char* filename, char mode, bool parallel);
  void file_close(hid_t file_id);
  void get_name(hid_t obj_id, char* name);
  int get_num_datasets(hid_t group_id);
  int get_num_groups(hid_t group_id);
  void get_datasets(hid_t group_id, char* name[]);
  void get_groups(hid_t group_id, char* name[]);
  void get_shape(hid_t obj_id, hsize_t* dims);
  void get_shape_attr(hid_t obj_id, const char* name, hsize_t* dims);
  bool object_exists(hid_t object_id, const char* name);
  hid_t open_dataset(hid_t group_id, const char* name);
  hid_t open_group(hid_t group_id, const char* name);
  void read_attr_double(hid_t obj_id, const char* name, double* buffer);
  void read_attr_int(hid_t obj_id, const char* name, int* buffer);
  void read_attr_string(hid_t obj_id, const char* name, size_t slen,
                        char* buffer);
  void read_complex(hid_t obj_id, const char* name,
                    std::complex<double>* buffer, bool indep);
  void read_double(hid_t obj_id, const char* name, double* buffer, bool indep);
  void read_int(hid_t obj_id, const char* name, int* buffer, bool indep);
  void read_llong(hid_t obj_id, const char* name, long long* buffer,
                  bool indep);
  void read_string(hid_t obj_id, const char* name, size_t slen, char* buffer,
                   bool indep);


  void read_tally_results(hid_t group_id, hsize_t n_filter, hsize_t n_score,
                          double* results);
  void write_attr_double(hid_t obj_id, int ndim, const hsize_t* dims,
                         const char* name, const double* buffer);
  void write_attr_int(hid_t obj_id, int ndim, const hsize_t* dims,
                      const char* name, const int* buffer);
  void write_attr_string(hid_t obj_id, const char* name, const char* buffer);
  void write_double(hid_t group_id, int ndim, const hsize_t* dims,
                    const char* name, const double* buffer, bool indep);
  void write_int(hid_t group_id, int ndim, const hsize_t* dims,
                 const char* name, const int* buffer, bool indep);
  void write_llong(hid_t group_id, int ndim, const hsize_t* dims,
                   const char* name, const long long* buffer, bool indep);
  void write_string(hid_t group_id, int ndim, const hsize_t* dims, size_t slen,
                    const char* name, char const* buffer, bool indep);
  void write_tally_results(hid_t group_id, hsize_t n_filter, hsize_t n_score,
                           const double* results);
} // extern "C"

//==============================================================================
// Template struct used to map types to HDF5 datatype IDs, which are stored
// using the type hid_t. By having a single static data member, the template can
// be specialized for each type we know of. The specializations appear in the
// .cpp file since they are definitions.
//==============================================================================

template<typename T>
struct H5TypeMap { static const hid_t type_id; };

//==============================================================================
// Templates/overloads for read_attribute
//==============================================================================

// Scalar version
template<typename T>
void read_attribute(hid_t obj_id, const char* name, T& buffer)
{
  read_attr(obj_id, name, H5TypeMap<T>::type_id, &buffer);
}

// array version
template<typename T, std::size_t N> inline void
read_attribute(hid_t obj_id, const char* name, std::array<T, N>& buffer)
{
  read_attr(obj_id, name, H5TypeMap<T>::type_id, buffer.data());
}

// vector version
template<typename T>
void read_attribute(hid_t obj_id, const char* name, std::vector<T>& vec)
{
  // Get shape of attribute array
  auto shape = attribute_shape(obj_id, name);

  // Allocate new array to read data into
  std::size_t size = 1;
  for (const auto x : shape)
    size *= x;
  vec.resize(size);

  // Read data from attribute
  read_attr(obj_id, name, H5TypeMap<T>::type_id, vec.data());
}

// Generic array version
template<typename T>
void read_attribute(hid_t obj_id, const char* name, xt::xarray<T>& arr)
{
  // Get shape of attribute array
  auto shape = attribute_shape(obj_id, name);

  // Allocate new array to read data into
  std::size_t size = 1;
  for (const auto x : shape)
    size *= x;
  std::vector<T> buffer(size);

  // Read data from attribute
  read_attr(obj_id, name, H5TypeMap<T>::type_id, buffer.data());

  // Adapt array into xarray
  arr = xt::adapt(buffer, shape);
}

// overload for std::string
inline void
read_attribute(hid_t obj_id, const char* name, std::string& str)
{
  // Create buffer to read data into
  auto n = attribute_typesize(obj_id, name);
  char* buffer = new char[n];

  // Read attribute and set string
  read_attr_string(obj_id, name, n, buffer);
  str = std::string{buffer, n};
  delete[] buffer;
}

// overload for std::vector<std::string>
inline void
read_attribute(hid_t obj_id, const char* name, std::vector<std::string>& vec)
{
  auto dims = attribute_shape(obj_id, name);
  auto m = dims[0];

  // Allocate a C char array to get strings
  auto n = attribute_typesize(obj_id, name);
  char* buffer = new char[m*n];

  // Read char data in attribute
  read_attr_string(obj_id, name, n, buffer);

  for (int i = 0; i < m; ++i) {
    // Determine proper length of string -- strlen doesn't work because
    // buffer[i] might not have any null characters
    std::size_t k = 0;
    for (; k < n; ++k) if (buffer[i*n + k] == '\0') break;

    // Create string based on (char*, size_t) constructor
    vec.emplace_back(&buffer[i*n], k);
  }
  delete[] buffer;
}

//==============================================================================
// Templates/overloads for read_dataset and related methods
//==============================================================================

// Template for scalars. We need to be careful that the compiler does not use
// this version of read_dataset for vectors, arrays, or other non-scalar types.
// enable_if_t allows us to conditionally remove the function from overload
// resolution when the type T doesn't meet a certain criterion.
template<typename T> inline
std::enable_if_t<std::is_scalar<std::decay_t<T>>::value>
read_dataset(hid_t obj_id, const char* name, T& buffer, bool indep=false)
{
  read_dataset_lowlevel(obj_id, name, H5TypeMap<T>::type_id, H5S_ALL, indep,
                        &buffer);
}

// overload for std::string
inline void
read_dataset(hid_t obj_id, const char* name, std::string& str, bool indep=false)
{
  // Create buffer to read data into
  auto n = dataset_typesize(obj_id, name);
  char* buffer = new char[n];

  // Read attribute and set string
  read_string(obj_id, name, n, buffer, indep);
  str = std::string{buffer, n};
}

// array version
template<typename T, std::size_t N> inline void
read_dataset(hid_t dset, const char* name, std::array<T, N>& buffer,
             bool indep=false)
{
  read_dataset_lowlevel(dset, name, H5TypeMap<T>::type_id, H5S_ALL, indep,
                        buffer.data());
}

// vector version
template <typename T>
void read_dataset(hid_t dset, std::vector<T>& vec, bool indep=false)
{
  // Get shape of dataset
  std::vector<hsize_t> shape = object_shape(dset);

  // Resize vector to appropriate size
  vec.resize(shape[0]);

  // Read data into vector
  read_dataset_lowlevel(dset, nullptr, H5TypeMap<T>::type_id, H5S_ALL, indep,
                        vec.data());
}

template <typename T>
void read_dataset(hid_t obj_id, const char* name, std::vector<T>& vec,
                  bool indep=false)
{
  hid_t dset = open_dataset(obj_id, name);
  read_dataset(dset, vec, indep);
  close_dataset(dset);
}

template <typename T>
void read_dataset(hid_t dset, xt::xarray<T>& arr, bool indep=false)
{
  // Get shape of dataset
  std::vector<hsize_t> shape = object_shape(dset);

  // Allocate space in the array to read data into
  std::size_t size = 1;
  for (const auto x : shape)
    size *= x;
  arr.resize(shape);

  // Read data from attribute
  read_dataset_lowlevel(dset, nullptr, H5TypeMap<T>::type_id, H5S_ALL, indep,
                        arr.data());
}

template<>
void read_dataset(hid_t dset, xt::xarray<std::complex<double>>& arr,
                  bool indep);

template <typename T>
void read_dataset(hid_t obj_id, const char* name, xt::xarray<T>& arr,
                  bool indep=false)
{
  // Open dataset and read array
  hid_t dset = open_dataset(obj_id, name);
  read_dataset(dset, arr, indep);
  close_dataset(dset);
}


template <typename T, std::size_t N>
void read_dataset(hid_t obj_id, const char* name, xt::xtensor<T, N>& arr,
                  bool indep=false)
{
  // Open dataset and read array
  hid_t dset = open_dataset(obj_id, name);

  // Get shape of dataset
  std::vector<hsize_t> hsize_t_shape = object_shape(dset);
  close_dataset(dset);

  // cast from hsize_t to size_t
  std::vector<size_t> shape(hsize_t_shape.size());
  for (int i = 0; i < shape.size(); i++) {
    shape[i] = static_cast<size_t>(hsize_t_shape[i]);
  }

  // Allocate new xarray to read data into
  xt::xarray<T> xarr(shape);

  // Read data from the dataset
  read_dataset(obj_id, name, xarr);

  // Copy into xtensor
  arr = xarr;
}

// overload for Position
inline void
read_dataset(hid_t obj_id, const char* name, Position& r, bool indep=false)
{
  std::array<double, 3> x;
  read_dataset(obj_id, name, x, indep);
  r.x = x[0];
  r.y = x[1];
  r.z = x[2];
}

template <typename T, std::size_t N>
inline void read_dataset_as_shape(hid_t obj_id, const char* name,
                                  xt::xtensor<T, N>& arr, bool indep=false)
{
  hid_t dset = open_dataset(obj_id, name);

  // Allocate new array to read data into
  std::size_t size = 1;
  for (const auto x : arr.shape())
    size *= x;
  std::vector<T> buffer(size);

  // Read data from attribute
  read_dataset_lowlevel(dset, nullptr, H5TypeMap<T>::type_id, H5S_ALL, indep,
                        buffer.data());

  // Adapt into xarray
  arr = xt::adapt(buffer, arr.shape());

  close_dataset(dset);
}


template <typename T, std::size_t N>
inline void read_nd_vector(hid_t obj_id, const char* name,
                           xt::xtensor<T, N>& result, bool must_have=false)
{
  if (object_exists(obj_id, name)) {
    read_dataset_as_shape(obj_id, name, result, true);
  } else if (must_have) {
    fatal_error(std::string("Must provide " + std::string(name) + "!"));
  }
}

//==============================================================================
// Templates/overloads for write_attribute
//==============================================================================

template<typename T> inline void
write_attribute(hid_t obj_id, const char* name, T buffer)
{
  write_attr(obj_id, 0, nullptr, name, H5TypeMap<T>::type_id, &buffer);
}

inline void
write_attribute(hid_t obj_id, const char* name, const char* buffer)
{
  write_attr_string(obj_id, name, buffer);
}

inline void
write_attribute(hid_t obj_id, const char* name, const std::string& buffer)
{
  write_attr_string(obj_id, name, buffer.c_str());
}

template<typename T, std::size_t N> inline void
write_attribute(hid_t obj_id, const char* name, const std::array<T, N>& buffer)
{
  hsize_t dims[] {N};
  write_attr(obj_id, 1, dims, name, H5TypeMap<T>::type_id, buffer.data());
}

template<typename T> inline void
write_attribute(hid_t obj_id, const char* name, const std::vector<T>& buffer)
{
  hsize_t dims[] {buffer.size()};
  write_attr(obj_id, 1, dims, name, H5TypeMap<T>::type_id, buffer.data());
}

inline void
write_attribute(hid_t obj_id, const char* name, Position r)
{
  std::array<double, 3> buffer {r.x, r.y, r.z};
  write_attribute(obj_id, name, buffer);
}



//==============================================================================
// Templates/overloads for write_dataset
//==============================================================================

// Template for scalars (ensured by SFINAE)
template<typename T> inline
std::enable_if_t<std::is_scalar<std::decay_t<T>>::value>
write_dataset(hid_t obj_id, const char* name, T buffer)
{
  write_dataset_lowlevel(obj_id, 0, nullptr, name, H5TypeMap<T>::type_id,
                         H5S_ALL, false, &buffer);
}

inline void
write_dataset(hid_t obj_id, const char* name, const char* buffer)
{
  write_string(obj_id, name, buffer, false);
}

template<typename T, std::size_t N> inline void
write_dataset(hid_t obj_id, const char* name, const std::array<T, N>& buffer)
{
  hsize_t dims[] {N};
  write_dataset_lowlevel(obj_id, 1, dims, name, H5TypeMap<T>::type_id,
                         H5S_ALL, false, buffer.data());
}

inline void
write_dataset(hid_t obj_id, const char* name,
              const std::vector<std::string>& buffer)
{
  auto n {buffer.size()};
  hsize_t dims[] {n};

  // Determine length of longest string, including \0
  size_t m = 1;
  for (const auto& s : buffer) {
    m = std::max(m, s.size() + 1);
  }

  // Copy data into contiguous buffer
  char* temp = new char[n*m];
  std::fill(temp, temp + n*m, '\0');
  for (int i = 0; i < n; ++i) {
    std::copy(buffer[i].begin(), buffer[i].end(), temp + i*m);
  }

  // Write 2D data
  write_string(obj_id, 1, dims, m, name, temp, false);

  // Free temp array
  delete[] temp;
}

template<typename T> inline void
write_dataset(hid_t obj_id, const char* name, const std::vector<T>& buffer)
{
  hsize_t dims[] {buffer.size()};
  write_dataset_lowlevel(obj_id, 1, dims, name, H5TypeMap<T>::type_id,
                         H5S_ALL, false, buffer.data());
}

// Template for xarray, xtensor, etc.
template<typename D> inline void
write_dataset(hid_t obj_id, const char* name, const xt::xcontainer<D>& arr)
{
  using T = typename D::value_type;
  auto s = arr.shape();
  std::vector<hsize_t> dims {s.cbegin(), s.cend()};
  write_dataset_lowlevel(obj_id, dims.size(), dims.data(), name,
                         H5TypeMap<T>::type_id, H5S_ALL, false, arr.data());
}

inline void
write_dataset(hid_t obj_id, const char* name, Position r)
{
  std::array<double, 3> buffer {r.x, r.y, r.z};
  write_dataset(obj_id, name, buffer);
}

inline void
write_dataset(hid_t obj_id, const char* name, std::string buffer)
{
  write_string(obj_id, name, buffer.c_str(), false);
}

} // namespace openmc
#endif // OPENMC_HDF5_INTERFACE_H
