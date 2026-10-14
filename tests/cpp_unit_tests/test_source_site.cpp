#include <cstddef>
#include <cstdint>
#include <cstring>

#include <catch2/catch_test_macros.hpp>

#include "openmc/particle_data.h"
#include "openmc/vector.h"

#ifdef OPENMC_MPI
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>

#include "openmc/initialize.h"
#include "openmc/message_passing.h"
#include "openmc/particle_type.h"
#endif

using namespace openmc;

#ifdef OPENMC_MPI

namespace {

// Set up MPI and the SourceSite datatype as openmc_init() does, once for the
// whole test run
class MPIListener : public Catch::EventListenerBase {
public:
  using Catch::EventListenerBase::EventListenerBase;

  void testRunStarting(const Catch::TestRunInfo&) override
  {
    initialize_mpi(MPI_COMM_WORLD);
  }

  void testRunEnded(const Catch::TestRunStats&) override
  {
    MPI_Type_free(&mpi::source_site);
    MPI_Type_free(&mpi::collision_track_site);
    MPI_Finalize();
  }
};

// Number of sites that each rank sends, and number of extra sites after them
// in each buffer, which must not be written to
constexpr int N_SITES {4};
constexpr int N_EXTRA {2};

// Site i of the given rank. Every byte of the site is first given a value
// that depends on the rank, the site and the position of the byte, so that
// members not set below also differ from site to site. Every member is then
// set to a value that differs from its default value and from the value of
// the same member in every other site.
SourceSite make_site(int rank, int i)
{
  unsigned char bytes[sizeof(SourceSite)];
  for (std::size_t k = 0; k < sizeof(SourceSite); ++k) {
    bytes[k] = static_cast<unsigned char>(1 + k + 37 * i + 101 * rank);
  }
  SourceSite site;
  std::memcpy(static_cast<void*>(&site), bytes, sizeof(SourceSite));

  int n = 100 * rank + 10 * i;
  double x = n;
  site.r = {x + 0.125, x + 0.25, x + 0.375};
  site.u = {x + 0.5, x + 0.625, x + 0.75};
  site.E = 1.0e6 + x;
  site.time = 1.0e-6 * (x + 1.0);
  site.wgt = x + 1.5;
  site.delayed_group = n + 1;
  site.surf_id = n + 2;
  site.particle = ParticleType {static_cast<int32_t>(n + 3)};
  site.parent_nuclide = n + 4;
  site.parent_id = (int64_t {1} << 40) + n + 5;
  site.progeny_id = (int64_t {1} << 41) + n + 6;
  site.wgt_born = x + 2.5;
  site.wgt_ww_born = x + 3.5;
  site.n_split = (int64_t {1} << 42) + n + 7;
  site.n_collision = n + 8;
  site.birth_mesh_bin = n + 9;
  return site;
}

// Fill every byte of a buffer of sites with the same value
void fill_bytes(vector<SourceSite>& sites, unsigned char value)
{
  std::memset(
    static_cast<void*>(sites.data()), value, sites.size() * sizeof(SourceSite));
}

} // namespace

CATCH_REGISTER_LISTENER(MPIListener)

TEST_CASE("The MPI datatype of SourceSite covers every byte of the struct")
{
  // Banks are sent as arrays of SourceSite, so the datatype must have the
  // extent of the struct
  MPI_Aint lb, extent;
  MPI_Type_get_extent(mpi::source_site, &lb, &extent);
  REQUIRE(lb == 0);
  REQUIRE(extent == static_cast<MPI_Aint>(sizeof(SourceSite)));

  // SourceSite has no padding: every byte belongs to a member. A datatype that
  // leaves out a member or part of one, or that gives a member a type of
  // another size, therefore does not hold as many bytes as the struct.
  MPI_Aint true_lb, true_extent;
  MPI_Type_get_true_extent(mpi::source_site, &true_lb, &true_extent);
  REQUIRE(true_lb == 0);
  REQUIRE(true_extent == static_cast<MPI_Aint>(sizeof(SourceSite)));
  int size;
  MPI_Type_size(mpi::source_site, &size);
  REQUIRE(size == static_cast<int>(sizeof(SourceSite)));
}

TEST_CASE("Source sites sent between ranks keep every member")
{
  // Each rank sends its sites to the next rank and receives the sites of the
  // previous one. With one rank, a rank sends its sites to itself.
  int next = (mpi::rank + 1) % mpi::n_procs;
  int previous = (mpi::rank + mpi::n_procs - 1) % mpi::n_procs;

  vector<SourceSite> sent(N_SITES + N_EXTRA);
  vector<SourceSite> received(N_SITES + N_EXTRA);
  fill_bytes(sent, 0xA5);
  fill_bytes(received, 0x5A);
  for (int i = 0; i < N_SITES; ++i) {
    sent[i] = make_site(mpi::rank, i);
  }

  MPI_Sendrecv(sent.data(), N_SITES, mpi::source_site, next, 0, received.data(),
    N_SITES, mpi::source_site, previous, 0, mpi::intracomm, MPI_STATUS_IGNORE);

  for (int i = 0; i < N_SITES; ++i) {
    const SourceSite expected = make_site(previous, i);
    const SourceSite& site = received[i];
    REQUIRE(site.r.x == expected.r.x);
    REQUIRE(site.r.y == expected.r.y);
    REQUIRE(site.r.z == expected.r.z);
    REQUIRE(site.u.x == expected.u.x);
    REQUIRE(site.u.y == expected.u.y);
    REQUIRE(site.u.z == expected.u.z);
    REQUIRE(site.E == expected.E);
    REQUIRE(site.time == expected.time);
    REQUIRE(site.wgt == expected.wgt);
    REQUIRE(site.delayed_group == expected.delayed_group);
    REQUIRE(site.surf_id == expected.surf_id);
    REQUIRE(site.particle.pdg_number() == expected.particle.pdg_number());
    REQUIRE(site.parent_nuclide == expected.parent_nuclide);
    REQUIRE(site.parent_id == expected.parent_id);
    REQUIRE(site.progeny_id == expected.progeny_id);
    REQUIRE(site.wgt_born == expected.wgt_born);
    REQUIRE(site.wgt_ww_born == expected.wgt_ww_born);
    REQUIRE(site.n_split == expected.n_split);
    REQUIRE(site.n_collision == expected.n_collision);
    REQUIRE(site.birth_mesh_bin == expected.birth_mesh_bin);

    // Since SourceSite has no padding, equal bytes mean equal members, also
    // for any member that is not compared above
    REQUIRE(std::memcmp(&site, &expected, sizeof(SourceSite)) == 0);
  }

  // Nothing was written past the received sites
  vector<SourceSite> untouched(N_EXTRA);
  fill_bytes(untouched, 0x5A);
  REQUIRE(std::memcmp(&received[N_SITES], untouched.data(),
            N_EXTRA * sizeof(SourceSite)) == 0);
}

#else

TEST_CASE("Source sites sent between ranks keep every member")
{
  SKIP("OpenMC was built without MPI");
}

#endif
