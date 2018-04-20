#include "initialize.h"

#include <cstddef>

#include "message_passing.h"
#include "openmc.h"


int openmc_init(const void* intracomm)
{
#ifdef OPENMC_MPI
  // Check if intracomm was passed
  MPI_Comm comm;
  if (intracomm) {
    comm = *static_cast<const MPI_Comm *>(intracomm);
  } else {
    comm = MPI_COMM_WORLD;
  }

  // Initialize MPI for C++
  openmc::initialize_mpi(comm);

  // Continue with rest of initialization
  MPI_Fint fcomm = MPI_Comm_c2f(comm);
  openmc_init_f(&fcomm);
#else
  openmc_init_f(nullptr);
#endif
  return 0;
}

namespace openmc {

#ifdef OPENMC_MPI
void initialize_mpi(MPI_Comm intracomm)
{
  openmc::mpi::intracomm = intracomm;

  // Initialize MPI
  int flag;
  if (!MPI_Initialized(&flag)) MPI_Init(nullptr, nullptr);

  // Determine number of processes and rank for each
  MPI_Comm_size(intracomm, &openmc::mpi::n_procs);
  MPI_Comm_rank(intracomm, &openmc::mpi::rank);

  // Set variable for Fortran side
  openmc_n_procs = openmc::mpi::n_procs;
  openmc_rank = openmc::mpi::rank;
  openmc_master = (openmc::mpi::rank == 0);

  // Create bank datatype
  Bank b;
  MPI_Aint disp[] {
    offsetof(Bank, wgt),
    offsetof(Bank, xyz),
    offsetof(Bank, uvw),
    offsetof(Bank, E),
    offsetof(Bank, delayed_group)
  };
  int blocks[] {1, 3, 3, 1, 1};
  MPI_Datatype types[] {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
  MPI_Type_create_struct(5, blocks, disp, types, &openmc::mpi::bank);
  MPI_Type_commit(&openmc::mpi::bank);
}

#endif // OPENMC_MPI
} // namespace openmc
