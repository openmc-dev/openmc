#include "openmc/message_passing.h"

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <vector>

using namespace openmc;

TEST_CASE("Double reductions span multiple MPI calls")
{
  int rank;
  int n_procs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
  REQUIRE(n_procs == 2);

  constexpr std::size_t count {8};
  std::vector<double> values(count);
  for (std::size_t i = 0; i < count; ++i) {
    values[i] = static_cast<double>(i + rank);
  }

  std::vector<double> reduced(rank == 0 ? count : 0);
  mpi::reduce(values.data(), rank == 0 ? reduced.data() : nullptr, count,
    MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    for (std::size_t i = 0; i < count; ++i) {
      REQUIRE(reduced[i] == 2.0 * i + 1.0);
    }
  }
}

TEST_CASE("In-place double reductions span multiple MPI calls")
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  constexpr std::size_t count {7};
  std::vector<double> values(count);
  for (std::size_t i = 0; i < count; ++i) {
    values[i] = static_cast<double>(i + rank);
  }

  if (rank == 0) {
    mpi::reduce(MPI_IN_PLACE, values.data(), count, MPI_SUM, 0, MPI_COMM_WORLD);
    for (std::size_t i = 0; i < count; ++i) {
      REQUIRE(values[i] == 2.0 * i + 1.0);
    }
  } else {
    mpi::reduce<double>(
      values.data(), nullptr, count, MPI_SUM, 0, MPI_COMM_WORLD);
  }
}

TEST_CASE("Integer reductions span multiple MPI calls")
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  constexpr std::size_t count {5};
  std::vector<int> values(count);
  for (std::size_t i = 0; i < count; ++i) {
    values[i] = static_cast<int>(i) + rank;
  }

  std::vector<int> reduced(rank == 0 ? count : 0);
  mpi::reduce(values.data(), rank == 0 ? reduced.data() : nullptr, count,
    MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    for (std::size_t i = 0; i < count; ++i) {
      REQUIRE(reduced[i] == 2 * static_cast<int>(i) + 1);
    }
  }
}

TEST_CASE("Reductions handle empty buffers")
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::vector<double> values;
  mpi::reduce(values.data(), rank == 0 ? values.data() : nullptr, 0, MPI_SUM, 0,
    MPI_COMM_WORLD);
}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int result = Catch::Session().run(argc, argv);
  MPI_Finalize();
  return result;
}
