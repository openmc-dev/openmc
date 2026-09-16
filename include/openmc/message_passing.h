#ifndef OPENMC_MESSAGE_PASSING_H
#define OPENMC_MESSAGE_PASSING_H

#include <cstddef>
#include <cstdint>

#ifdef OPENMC_MPI
#include <mpi.h>
#endif

#include "openmc/vector.h"

namespace openmc {
namespace mpi {

extern int rank;
extern int n_procs;
extern bool master;

#ifdef OPENMC_MPI
extern MPI_Datatype source_site;
extern MPI_Datatype collision_track_site;
extern MPI_Comm intracomm;

//==============================================================================
// Template struct used to map types to MPI datatypes
// By having a single static data member, the template can
// be specialized for each type we know of. The specializations appear in the
// .cpp file since they are definitions.
//==============================================================================
template<typename T>
struct MPITypeMap {
  static const MPI_Datatype mpi_type;
};

//! Reduce an array without narrowing the element count.
//
//! \param sendbuf Input buffer, or MPI_IN_PLACE on the root process
//! \param recvbuf Output buffer on the root process; ignored on other processes
//! \param count Number of elements in each buffer
//! \param op Operation to apply during the reduction
//! \param root Rank of the root process
//! \param comm Communicator over which to perform the reduction
void reduce_buffer(const void* sendbuf, void* recvbuf, std::size_t count,
  MPI_Datatype datatype, std::size_t type_size, MPI_Op op, int root,
  MPI_Comm comm);

template<typename T>
void reduce(const void* sendbuf, T* recvbuf, std::size_t count, MPI_Op op,
  int root, MPI_Comm comm)
{
  reduce_buffer(sendbuf, recvbuf, count, MPITypeMap<T>::mpi_type, sizeof(T), op,
    root, comm);
}
#endif

// Calculates global indices of the bank particles
// across all ranks using a parallel scan. This is used to write
// the surface source file in parallel runs. It will probably
// be used in the future for other types of bank like particles
// in flight used to kick off transient simulations.
//
// More abstractly, this just takes a number from each MPI rank,
// and returns a vector which is the exclusive parallel scan across
// all of those numbers, having a length of the number of MPI ranks
// plus one.
vector<int64_t> calculate_parallel_index_vector(int64_t size);

} // namespace mpi
} // namespace openmc

#endif // OPENMC_MESSAGE_PASSING_H
