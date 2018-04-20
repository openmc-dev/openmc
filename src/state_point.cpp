#include "state_point.h"

#include <algorithm>
#include <vector>

#include "mpi.h"
#include "message_passing.h"
#include "openmc.h"

namespace openmc {

void
write_source_bank(hid_t group_id, int64_t* work_index, const Bank* source_bank)
{
  // Create type for array of 3 reals
  hsize_t dims[] {3};
  hid_t triplet = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dims);

  // Create bank datatype
  hid_t banktype = H5Tcreate(H5T_COMPOUND, sizeof(struct Bank));
  H5Tinsert(banktype, "wgt", HOFFSET(Bank, wgt), H5T_NATIVE_DOUBLE);
  H5Tinsert(banktype, "xyz", HOFFSET(Bank, xyz), triplet);
  H5Tinsert(banktype, "uvw", HOFFSET(Bank, uvw), triplet);
  H5Tinsert(banktype, "E", HOFFSET(Bank, E), H5T_NATIVE_DOUBLE);
  H5Tinsert(banktype, "delayed_group", HOFFSET(Bank, delayed_group), H5T_NATIVE_INT);

#ifdef PHDF5
  // Set size of total dataspace for all procs and rank
  dims[0] = n_particles;
  hid_t dspace = H5Screate_simple(1, dims, H5P_DEFAULT);
  hid_t dset = H5Dcreate(group_id, "source_bank", banktype, dspace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Create another data space but for each proc individually
  hsize_t count[] {openmc_work};
  hid_t memspace = H5Screate_simple(1, count, H5P_DEFAULT);

  // Select hyperslab for this dataspace
  hsize_t start[] {work_index[openmc::mpi::rank]};
  H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, nullptr, count, nullptr);

  // Set up the property list for parallel writing
  hid_t plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE);

  // Write data to file in parallel
  H5Dwrite(dset, banktype, memspace, dspace, memspace, plist, source_bank);

  // Free resources
  H5Sclose(dspace);
  h5sclose(memspace);
  H5Dclose(dset);
  H5Pclose(plist);

#else

  if (openmc_master) {
    // Create dataset big enough to hold all source sites
    hsize_t dims[] {static_cast<hsize_t>(n_particles)};
    hid_t dspace = H5Screate_simple(1, dims, H5P_DEFAULT);
    hid_t dset = H5Dcreate(group_id, "source_bank", banktype, dspace,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Save source bank sites since the souce_bank array is overwritten below
#ifdef OPENMC_MPI
    std::vector<double> temp_source {source_bank, source_bank + openmc_work};
#endif

    for (int i = 0; i < openmc::mpi::n_procs; ++i) {
      // Create memory space
      hsize_t count[] {static_cast<hsize_t>(work_index[i+1] - work_index[i])};
      hid_t memspace = H5Screate_simple(1, count, H5P_DEFAULT);

#ifdef OPENMC_MPI
      // Receive source sites from other processes
      if (i > 0)
        MPI_Recv(source_bank, count[0], openmc::mpi::bank, i, i,
                 openmc::mpi::intracomm, MPI_STATUS_IGNORE);
#endif

      // Select hyperslab for this dataspace
      dspace = H5Dget_space(dset);
      hsize_t start[] {static_cast<hsize_t>(work_index[i])};
      H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, nullptr, count, nullptr);

      // Write data to hyperslab
      H5Dwrite(dset, banktype, memspace, dspace, H5P_DEFAULT, source_bank);

      H5Sclose(memspace);
      H5Sclose(dspace);
    }

    // Close all ids
    H5Dclose(dset);

#ifdef OPENMC_MPI
    // Restore state of source bank
    std::copy(temp_source.begin(), temp_source.end(), source_bank);
#endif
  } else {
#ifdef OPENMC_MPI
    MPI_Send(source_bank, openmc_work, openmc::mpi::bank, 0, openmc::mpi::rank,
             openmc::mpi::mpi_intracomm);
#endif
  }
#endif

  H5Tclose(banktype);
  H5Tclose(triplet);
}

} // namespace openmc
