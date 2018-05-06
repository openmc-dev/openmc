#ifndef HDF5_INTERFACE_H
#define HDF5_INTERFACE_H

#include "hdf5.h"
#include "hdf5_hl.h"

#include <array>
#include <string>
#include <sstream>
#include <complex.h>


namespace openmc {

extern "C" bool attribute_exists(hid_t obj_id, const char* name);
extern "C" size_t attribute_typesize(hid_t obj_id, const char* name);
extern "C" hid_t create_group(hid_t parent_id, const char* name);
hid_t create_group(hid_t parent_id, const std::string& name);
extern "C" void close_dataset(hid_t dataset_id);
extern "C" void close_group(hid_t group_id);
extern "C" int dataset_ndims(hid_t dset);
extern "C" size_t dataset_typesize(hid_t dset);
extern "C" hid_t file_open(const char* filename, char mode, bool parallel);
hid_t file_open(const std::string& filename, char mode, bool parallel);
extern "C" void file_close(hid_t file_id);
extern "C" void get_name(hid_t obj_id, char* name);
extern "C" int get_num_datasets(hid_t group_id);
extern "C" int get_num_groups(hid_t group_id);
extern "C" void get_datasets(hid_t group_id, char* name[]);
extern "C" void get_groups(hid_t group_id, char* name[]);
extern "C" void get_shape(hid_t obj_id, hsize_t* dims);
extern "C" void get_shape_attr(hid_t obj_id, const char* name, hsize_t* dims);
extern "C" bool object_exists(hid_t object_id, const char* name);
extern "C" hid_t open_dataset(hid_t group_id, const char* name);
extern "C" hid_t open_group(hid_t group_id, const char* name);
bool using_mpio_device(hid_t obj_id);

void read_attr(hid_t obj_id, const char* name, hid_t mem_type_id,
               const void* buffer);
extern "C" void read_attr_double(hid_t obj_id, const char* name, double* buffer);
extern "C" void read_attr_int(hid_t obj_id, const char* name, int* buffer);
extern "C" void read_attr_string(hid_t obj_id, const char* name, size_t slen,
                                 char* buffer);

void read_dataset(hid_t obj_id, const char* name, hid_t mem_type_id,
                  void* buffer, bool indep);
extern "C" void read_double(hid_t obj_id, const char* name, double* buffer,
                            bool indep);
extern "C" void read_int(hid_t obj_id, const char* name, int* buffer,
                         bool indep);
extern "C" void read_llong(hid_t obj_id, const char* name, long long* buffer,
                           bool indep);
extern "C" void read_string(hid_t obj_id, const char* name, size_t slen,
                            char* buffer, bool indep);
extern "C" void read_complex(hid_t obj_id, const char* name,
                             double _Complex* buffer, bool indep);

extern "C" void read_tally_results(hid_t group_id, hsize_t n_filter,
                                   hsize_t n_score, double* results);

void write_attr(hid_t obj_id, int ndim, const hsize_t* dims, const char* name,
                hid_t mem_type_id, const void* buffer);
extern "C" void write_attr_double(hid_t obj_id, int ndim, const hsize_t* dims,
                                  const char* name, const double* buffer);
extern "C" void write_attr_int(hid_t obj_id, int ndim, const hsize_t* dims,
                               const char* name, const int* buffer);
extern "C" void write_attr_string(hid_t obj_id, const char* name, const char* buffer);


void write_dataset(hid_t group_id, int ndim, const hsize_t* dims, const char* name,
                   hid_t mem_type_id, const void* buffer, bool indep);
extern "C" void write_double(hid_t group_id, int ndim, const hsize_t* dims,
                             const char* name, const double* buffer, bool indep);
extern "C" void write_int(hid_t group_id, int ndim, const hsize_t* dims,
                          const char* name, const int* buffer, bool indep);
extern "C" void write_llong(hid_t group_id, int ndim, const hsize_t* dims,
                            const char* name, const long long* buffer, bool indep);

extern "C" void write_string(hid_t group_id, int ndim, const hsize_t* dims, size_t slen,
                             const char* name, char const* buffer, bool indep);
void write_string(hid_t group_id, const char* name, const std::string& buffer, bool indep);


extern "C" void write_tally_results(hid_t group_id, hsize_t n_filter, hsize_t n_score,
                                    const double* results);

template<std::size_t array_len> void
write_double_1D(hid_t group_id, char const *name,
                std::array<double, array_len> &buffer)
{
  hsize_t dims[1] {array_len};
  write_dataset(group_id, 1, dims, name, H5T_NATIVE_DOUBLE,
                buffer.data(), false);
}

} // namespace openmc
#endif //HDF5_INTERFACE_H
