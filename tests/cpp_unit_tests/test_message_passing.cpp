#include <climits>
#include <cstddef>
#include <cstdint>

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>

#include "openmc/message_passing.h"
#include "openmc/vector.h"

using namespace openmc;

template<typename T>
void check_broadcast(int root)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  for (std::size_t count : {0, 1, 4, 9}) {
    vector<T> buffer(count, T {-1});
    if (rank == root) {
      for (std::size_t i = 0; i < count; ++i) {
        buffer[i] = static_cast<T>(i + 1);
      }
    }
    mpi::broadcast(buffer.data(), buffer.size(), root, MPI_COMM_WORLD);
    for (std::size_t i = 0; i < count; ++i) {
      CHECK(buffer[i] == static_cast<T>(i + 1));
    }
  }
}

TEST_CASE("Broadcast contiguous buffers from either rank")
{
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  for (int root : {0, size - 1}) {
    check_broadcast<int>(root);
    check_broadcast<int64_t>(root);
    check_broadcast<double>(root);
  }
}

TEST_CASE("Broadcast counts exceeding the legacy MPI limit")
{
  // A zero-size datatype exercises large counts without allocating a large
  // buffer. The legacy path must split the count rather than narrow it.
  MPI_Datatype empty;
  REQUIRE(MPI_Type_contiguous(0, MPI_BYTE, &empty) == MPI_SUCCESS);
  REQUIRE(MPI_Type_commit(&empty) == MPI_SUCCESS);
  int sentinel = 42;
  for (std::size_t count :
    {static_cast<std::size_t>(INT_MAX), static_cast<std::size_t>(INT_MAX) + 1,
      2 * static_cast<std::size_t>(INT_MAX) + 1}) {
    mpi::broadcast_buffer(&sentinel, count, empty, 0, 0, MPI_COMM_WORLD);
    CHECK(sentinel == 42);
  }
  REQUIRE(MPI_Type_free(&empty) == MPI_SUCCESS);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int result = Catch::Session().run(argc, argv);
  MPI_Finalize();
  return result;
}
