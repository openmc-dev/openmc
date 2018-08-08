#ifndef OPENMC_HDF5_INTERFACE_H
#define OPENMC_HDF5_INTERFACE_H

#include "hdf5.h"
#include "hdf5_hl.h"

#include <array>
#include <cstddef>
#include <string>
#include <sstream>
#include <vector>
#include <complex.h>

#include "position.h"


namespace openmc {

//==============================================================================
// Low-level internal functions
//==============================================================================

void read_attr(hid_t obj_id, const char* name, hid_t mem_type_id,
               const void* buffer);
void write_attr(hid_t obj_id, int ndim, const hsize_t* dims, const char* name,
                hid_t mem_type_id, const void* buffer);
void read_dataset(hid_t obj_id, const char* name, hid_t mem_type_id,
                  void* buffer, bool indep);
void write_dataset(hid_t group_id, int ndim, const hsize_t* dims, const char* name,
                   hid_t mem_type_id, const void* buffer, bool indep);
bool using_mpio_device(hid_t obj_id);

//==============================================================================
// Normal functions that are used to read/write files
//==============================================================================

hid_t create_group(hid_t parent_id, const std::string& name);
hid_t file_open(const std::string& filename, char mode, bool parallel=false);
void write_string(hid_t group_id, const char* name, const std::string& buffer,
                  bool indep);

void
read_nd_vector(hid_t obj_id, const char* name, std::vector<double>& result,
               bool must_have = false);

void
read_nd_vector(hid_t obj_id, const char* name,
               std::vector<std::vector<double> >& result,
               bool must_have = false);

void
read_nd_vector(hid_t obj_id, const char* name,
               std::vector<std::vector<int> >& result, bool must_have = false);

void
read_nd_vector(hid_t obj_id, const char* name,
               std::vector<std::vector<std::vector<double> > >& result,
               bool must_have = false);

void
read_nd_vector(hid_t obj_id, const char* name,
               std::vector<std::vector<std::vector<int> > >& result,
               bool must_have = false);

void
read_nd_vector(hid_t obj_id, const char* name,
               std::vector<std::vector<std::vector<std::vector<double> > > >& result,
               bool must_have = false);

void
read_nd_vector(hid_t obj_id, const char* name,
               std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > >& result,
               bool must_have = false);

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
  size_t dataset_typesize(hid_t dset);
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
                             double _Complex* buffer, bool indep);
  void read_double(hid_t obj_id, const char* name, double* buffer,
                            bool indep);
  void read_int(hid_t obj_id, const char* name, int* buffer,
                         bool indep);
  void read_llong(hid_t obj_id, const char* name, long long* buffer,
                           bool indep);
  void read_string(hid_t obj_id, const char* name, size_t slen,
                            char* buffer, bool indep);


  void read_tally_results(hid_t group_id, hsize_t n_filter,
                                   hsize_t n_score, double* results);
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
// Template functions used to provide simple interface to lower-level functions
//==============================================================================

template<typename T> inline void
write_attribute(hid_t obj_id, const char* name, T buffer)
{
  write_attr(obj_id, name, 0, nullptr, H5TypeMap<T>::type_id, &buffer);
}

template<> inline void
write_attribute<const char*>(hid_t obj_id, const char* name, const char* buffer)
{
  write_attr_string(obj_id, name, buffer);
}

template<typename T, std::size_t N> inline void
write_attribute(hid_t obj_id, const char* name, const std::array<T, N>& buffer)
{
  hsize_t dims[] {N};
  write_attr(obj_id, 1, dims, name, H5TypeMap<T>::type_id, buffer.data());
}

template<typename T> inline void
write_dataset(hid_t obj_id, const char* name, T buffer)
{
  write_dataset(obj_id, 0, nullptr, name, H5TypeMap<T>::type_id, &buffer, false);
}

template<> inline void
write_dataset<const char*>(hid_t obj_id, const char* name, const char* buffer)
{
  write_string(obj_id, name, buffer, false);
}

template<typename T, std::size_t N> inline void
write_dataset(hid_t obj_id, const char* name, const std::array<T, N>& buffer)
{
  hsize_t dims[] {N};
  write_dataset(obj_id, 1, dims, name, H5TypeMap<T>::type_id, buffer.data(), false);
}

inline void
write_dataset(hid_t obj_id, const char* name, Position r)
{
  std::array<double, 3> buffer {r.x, r.y, r.z};
  write_dataset(obj_id, name, buffer);
}

} // namespace openmc
#endif // OPENMC_HDF5_INTERFACE_H
