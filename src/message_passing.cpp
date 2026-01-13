#include "openmc/message_passing.h"

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

} // namespace mpi

} // namespace openmc
