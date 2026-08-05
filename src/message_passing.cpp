#include "openmc/message_passing.h"

#include <algorithm>
#include <climits>
#include <limits>

namespace openmc {
namespace mpi {

int rank {0};
int n_procs {1};
bool master {true};

#ifdef OPENMC_MPI
MPI_Comm intracomm {MPI_COMM_NULL};
MPI_Datatype source_site {MPI_DATATYPE_NULL};
MPI_Datatype collision_track_site {MPI_DATATYPE_NULL};
#endif

extern "C" bool openmc_master()
{
  return mpi::master;
}

vector<int64_t> calculate_parallel_index_vector(int64_t size)
{
  vector<int64_t> result;
  result.resize(n_procs + 1);
  result[0] = 0;

#ifdef OPENMC_MPI

  // Populate the result with cumulative sum of the number of
  // surface source banks per process
  int64_t scan_total;
  MPI_Scan(&size, &scan_total, 1, MPI_INT64_T, MPI_SUM, intracomm);
  MPI_Allgather(
    &scan_total, 1, MPI_INT64_T, result.data() + 1, 1, MPI_INT64_T, intracomm);
#else
  result[1] = size;
#endif

  return result;
}

#ifdef OPENMC_MPI
void reduce_buffer(const void* sendbuf, void* recvbuf, std::size_t count,
  MPI_Datatype datatype, std::size_t type_size, MPI_Op op, int root,
  MPI_Comm comm)
{
  // Determine the maximum number of elements accepted by the selected API.
#ifdef OPENMC_HAVE_MPI_REDUCE_C
  const std::size_t chunk_limit {
    static_cast<std::size_t>(std::numeric_limits<MPI_Count>::max())};
#else
  constexpr std::size_t chunk_limit {INT_MAX};
#endif

  // A legacy MPI reduction accepts an int count, so split the logical buffer
  // into bounded chunks before converting the count.
  for (std::size_t offset = 0; offset < count; offset += chunk_limit) {
    const std::size_t chunk_size = std::min(count - offset, chunk_limit);

    // MPI_IN_PLACE must be passed unchanged on the root. Other buffers are
    // advanced by the chunk offset in bytes because their element type is
    // represented by an MPI datatype rather than a C++ pointer type here.
    const void* send_chunk =
      sendbuf == MPI_IN_PLACE
        ? MPI_IN_PLACE
        : static_cast<const char*>(sendbuf) + offset * type_size;
    void* recv_chunk = recvbuf == nullptr
                         ? nullptr
                         : static_cast<char*>(recvbuf) + offset * type_size;

#ifdef OPENMC_HAVE_MPI_REDUCE_C
    MPI_Reduce_c(send_chunk, recv_chunk, static_cast<MPI_Count>(chunk_size),
      datatype, op, root, comm);
#else
    MPI_Reduce(send_chunk, recv_chunk, static_cast<int>(chunk_size), datatype,
      op, root, comm);
#endif
  }
}

// Specializations of the MPITypeMap template struct
template<>
const MPI_Datatype MPITypeMap<int>::mpi_type = MPI_INT;
template<>
const MPI_Datatype MPITypeMap<double>::mpi_type = MPI_DOUBLE;
template<>
const MPI_Datatype MPITypeMap<int64_t>::mpi_type = MPI_INT64_T;
#endif

} // namespace mpi

} // namespace openmc
