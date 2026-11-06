#include <cstdint>

#include <catch2/catch_test_macros.hpp>

#include "openmc/bank.h"
#include "openmc/constants.h"
#include "openmc/geometry.h"
#include "openmc/message_passing.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/particle_type.h"
#include "openmc/settings.h"
#include "openmc/shared_array.h"

#ifdef OPENMC_MPI
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>

#include "openmc/initialize.h"
#endif

using namespace openmc;

namespace {

#ifdef OPENMC_MPI
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
#endif

// Settings of a fixed-source run with a shared secondary bank. Particle
// positions are stored in the coordinate levels, so each particle gets one.
class SharedSecondaryBank {
public:
  SharedSecondaryBank()
    : n_coord_levels_ {model::n_coord_levels}, run_mode_ {settings::run_mode},
      shared_ {settings::use_shared_secondary_bank}
  {
    model::n_coord_levels = 1;
    settings::run_mode = RunMode::FIXED_SOURCE;
    settings::use_shared_secondary_bank = true;
  }

  ~SharedSecondaryBank()
  {
    model::n_coord_levels = n_coord_levels_;
    settings::run_mode = run_mode_;
    settings::use_shared_secondary_bank = shared_;
  }

private:
  int n_coord_levels_;
  RunMode run_mode_;
  bool shared_;
};

// Each state of the parent neutron banks these secondaries, in this order
constexpr int N_PARENT_STATES {4};
const ParticleType SECONDARY_TYPES[] {ParticleType::neutron(),
  ParticleType::neutron(), ParticleType::photon(), ParticleType::electron(),
  ParticleType::positron()};
constexpr int N_SECONDARIES {5};

// Clock of the parent neutron in state i. The lifetime differs from the time,
// which also counts the time before the birth of the neutron's generation.
double parent_lifetime(int64_t i)
{
  return 1.25e-6 * (i + 1);
}

double parent_time(int64_t i)
{
  return 3.0e-6 + parent_lifetime(i);
}

} // namespace

#ifdef OPENMC_MPI
CATCH_REGISTER_LISTENER(MPIListener)

TEST_CASE("The MPI datatype of SourceSite spans the whole struct")
{
  // Banks are sent as arrays of SourceSite, so the datatype must have the
  // extent of the struct
  MPI_Aint lb, extent;
  MPI_Type_get_extent(mpi::source_site, &lb, &extent);
  REQUIRE(lb == 0);
  REQUIRE(extent == static_cast<MPI_Aint>(sizeof(SourceSite)));
}
#endif

TEST_CASE_METHOD(SharedSecondaryBank,
  "Neutrons revived from the shared secondary bank keep their lifetime")
{
  // All secondaries are banked on the first rank, so with more than one rank
  // the bank is redistributed across ranks before they are transported
  SharedArray<SourceSite> bank;
  if (mpi::rank == 0) {
    Particle p;
    p.type() = ParticleType::neutron();
    p.wgt() = 1.0;
    p.r() = {1.0, 2.0, 3.0};
    p.u() = {0.0, 0.0, 1.0};
    p.E() = 2.0e6;
    for (int i = 0; i < N_PARENT_STATES; ++i) {
      p.time() = parent_time(i);
      p.lifetime() = parent_lifetime(i);
      // A neutron split by a weight window
      p.split(0.5);
      // An extra neutron of an (n,2n) reaction, then the other particles
      for (int k = 1; k < N_SECONDARIES; ++k) {
        CHECK(p.create_secondary(1.0, p.u(), 1.0e6, SECONDARY_TYPES[k]));
      }
    }
    for (const auto& site : p.local_secondary_bank()) {
      bank.thread_unsafe_append(site);
    }
  }

  const int64_t n_sites = N_PARENT_STATES * N_SECONDARIES;
  int64_t n_total = synchronize_global_secondary_bank(bank);
  REQUIRE(n_total == n_sites);

  // Every rank receives its share of the bank
  int64_t share =
    n_sites / mpi::n_procs + (mpi::rank < n_sites % mpi::n_procs ? 1 : 0);
  REQUIRE(bank.size() == share);

  // One particle revives all sites, as in the transport loop
  Particle q;
  q.lifetime() = -1.0;
  for (int64_t j = 0; j < bank.size(); ++j) {
    const auto& site = bank[j];
    int64_t i = site.progeny_id / N_SECONDARIES;
    int k = site.progeny_id % N_SECONDARIES;
    REQUIRE(site.progeny_id >= 0);
    REQUIRE(i < N_PARENT_STATES);

    q.event_revive_from_secondary(site);
    REQUIRE(q.type() == SECONDARY_TYPES[k]);
    REQUIRE(q.time() == parent_time(i));
    if (q.type().is_neutron()) {
      REQUIRE(q.lifetime() == parent_lifetime(i));
    } else {
      REQUIRE(q.lifetime() == 0.0);
    }
  }

  // A site that is not a secondary, such as a fission site or an external
  // source site, starts a new clock
  SourceSite site;
  site.particle = ParticleType::neutron();
  site.r = {1.0, 2.0, 3.0};
  site.u = {0.0, 0.0, 1.0};
  site.E = 2.0e6;
  site.time = 5.0e-6;
  q.lifetime() = 1.0e-6;
  q.event_revive_from_secondary(site);
  REQUIRE(q.time() == 5.0e-6);
  REQUIRE(q.lifetime() == 0.0);
}

TEST_CASE_METHOD(SharedSecondaryBank,
  "Neutrons revived from the shared secondary bank keep their delayed group")
{
  // Neutrons born from precursors of delayed groups 1 to N_PARENT_STATES bank
  // the secondaries of SECONDARY_TYPES on the first rank. Then a photon banks a
  // split photon and a neutron. It is given a nonzero delayed group, which
  // transport never gives a photon, so that a copy of it would show.
  const int64_t n_neutron_parent_sites = N_PARENT_STATES * N_SECONDARIES;
  const int64_t n_sites = n_neutron_parent_sites + 2;
  SharedArray<SourceSite> bank;
  if (mpi::rank == 0) {
    Particle p;
    p.type() = ParticleType::neutron();
    p.wgt() = 1.0;
    p.r() = {1.0, 2.0, 3.0};
    p.u() = {0.0, 0.0, 1.0};
    p.E() = 2.0e6;
    for (int i = 0; i < N_PARENT_STATES; ++i) {
      p.delayed_group() = i + 1;
      p.split(0.5);
      for (int k = 1; k < N_SECONDARIES; ++k) {
        CHECK(p.create_secondary(1.0, p.u(), 1.0e6, SECONDARY_TYPES[k]));
      }
    }
    p.type() = ParticleType::photon();
    p.delayed_group() = N_PARENT_STATES + 1;
    p.split(0.5);
    CHECK(p.create_secondary(1.0, p.u(), 1.0e6, ParticleType::neutron()));
    for (const auto& site : p.local_secondary_bank()) {
      bank.thread_unsafe_append(site);
    }
  }

  int64_t n_total = synchronize_global_secondary_bank(bank);
  REQUIRE(n_total == n_sites);
  int64_t share =
    n_sites / mpi::n_procs + (mpi::rank < n_sites % mpi::n_procs ? 1 : 0);
  REQUIRE(bank.size() == share);

  Particle q;
  for (int64_t j = 0; j < bank.size(); ++j) {
    const auto& site = bank[j];
    int64_t id = site.progeny_id;
    REQUIRE(id >= 0);
    REQUIRE(id < n_sites);

    // Only the neutrons banked by a neutron take its delayed group
    ParticleType type;
    int group = 0;
    if (id < n_neutron_parent_sites) {
      type = SECONDARY_TYPES[id % N_SECONDARIES];
      if (type.is_neutron()) {
        group = static_cast<int>(id / N_SECONDARIES) + 1;
      }
    } else if (id == n_neutron_parent_sites) {
      type = ParticleType::photon();
    } else {
      type = ParticleType::neutron();
    }
    REQUIRE(site.particle == type);
    REQUIRE(site.delayed_group == group);

    q.delayed_group() = -1;
    q.event_revive_from_secondary(site);
    REQUIRE(q.type() == type);
    REQUIRE(q.delayed_group() == group);
  }

  // A site that is not a secondary, such as a fission site, sets its own
  // delayed group
  SourceSite site;
  site.particle = ParticleType::neutron();
  site.r = {1.0, 2.0, 3.0};
  site.u = {0.0, 0.0, 1.0};
  site.E = 2.0e6;
  site.delayed_group = 2;
  q.delayed_group() = -1;
  q.event_revive_from_secondary(site);
  REQUIRE(q.delayed_group() == 2);
}
