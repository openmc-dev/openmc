#include "openmc/message_passing.h"

namespace openmc {
namespace mpi {

int rank {0};
int n_procs {1};
bool master {true};

#ifdef OPENMC_MPI
MPI_Comm intracomm {MPI_COMM_NULL};
MPI_Datatype source_site {MPI_DATATYPE_NULL};
#endif

extern "C" bool openmc_master()
{
  return mpi::master;
}

vector<int64_t> calculate_parallel_index_vector(const int64_t size)
{
  vector<int64_t> result;
  result.reserve(n_procs + 1);

#ifdef OPENMC_MPI
  result.resize(n_procs);
  result[0] = 0;
  vector<int64_t> bank_size(n_procs + 1);

  // Populate the result with cumulative sum of the number of
  // surface source banks per process
  MPI_Scan(&size, bank_size.data(), 1, MPI_INT64_T, MPI_SUM, intracomm);
  MPI_Allgather(bank_size.data() + 1, 1, MPI_INT64_T, result.data(), 1,
    MPI_INT64_T, intracomm);
#else
  result.push_back(0);
  result.push_back(size);
#endif

  return result;
}

} // namespace mpi

} // namespace openmc
