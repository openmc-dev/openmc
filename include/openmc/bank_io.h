#ifndef OPENMC_BANK_IO_H
#define OPENMC_BANK_IO_H

#include "hdf5.h"

#include "openmc/message_passing.h"
#include "openmc/span.h"
#include "openmc/vector.h"

#include <algorithm>

#ifdef OPENMC_MPI
#include <mpi.h>
#endif

namespace openmc {

template <typename SiteType>
void write_bank_dataset(const char* dataset_name, hid_t group_id,
  span<SiteType> bank, const vector<int64_t>& bank_index, hid_t banktype
#ifdef OPENMC_MPI
  , MPI_Datatype mpi_dtype
#endif
)
{
  int64_t dims_size = bank_index.back();
  int64_t count_size = bank_index[mpi::rank + 1] - bank_index[mpi::rank];

#ifdef PHDF5
  hsize_t dims[] {static_cast<hsize_t>(dims_size)};
  hid_t dspace = H5Screate_simple(1, dims, nullptr);
  hid_t dset = H5Dcreate(
    group_id, dataset_name, banktype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hsize_t count[] {static_cast<hsize_t>(count_size)};
  hid_t memspace = H5Screate_simple(1, count, nullptr);

  hsize_t start[] {static_cast<hsize_t>(bank_index[mpi::rank])};
  H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, nullptr, count, nullptr);

  hid_t plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

  H5Dwrite(dset, banktype, memspace, dspace, plist, bank.data());

  H5Sclose(dspace);
  H5Sclose(memspace);
  H5Dclose(dset);
  H5Pclose(plist);
#else
  if (mpi::master) {
    hsize_t dims[] {static_cast<hsize_t>(dims_size)};
    hid_t dspace = H5Screate_simple(1, dims, nullptr);
    hid_t dset = H5Dcreate(
      group_id, dataset_name, banktype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

#ifdef OPENMC_MPI
    vector<SiteType> temp_bank {bank.begin(), bank.end()};
#endif

    for (int i = 0; i < mpi::n_procs; ++i) {
      hsize_t count[] {static_cast<hsize_t>(bank_index[i + 1] - bank_index[i])};
      hid_t memspace = H5Screate_simple(1, count, nullptr);

#ifdef OPENMC_MPI
      if (i > 0) {
        MPI_Recv(bank.data(), count[0], mpi_dtype, i, i, mpi::intracomm,
          MPI_STATUS_IGNORE);
      }
#endif

      hid_t dspace_rank = H5Dget_space(dset);
      hsize_t start[] {static_cast<hsize_t>(bank_index[i])};
      H5Sselect_hyperslab(
        dspace_rank, H5S_SELECT_SET, start, nullptr, count, nullptr);

      H5Dwrite(dset, banktype, memspace, dspace_rank, H5P_DEFAULT, bank.data());

      H5Sclose(memspace);
      H5Sclose(dspace_rank);
    }

    H5Dclose(dset);

#ifdef OPENMC_MPI
    std::copy(temp_bank.begin(), temp_bank.end(), bank.begin());
#endif
  }
#ifdef OPENMC_MPI
  else {
    if (!bank.empty()) {
      MPI_Send(
        bank.data(), bank.size(), mpi_dtype, 0, mpi::rank, mpi::intracomm);
    }
  }
#endif
#endif
}

} // namespace openmc

#endif // OPENMC_BANK_IO_H
