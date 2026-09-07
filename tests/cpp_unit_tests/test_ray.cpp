#include <algorithm>

#include <catch2/catch_test_macros.hpp>

#include "openmc/geometry.h"
#include "openmc/particle.h"
#include "openmc/particle_type.h"
#include "openmc/random_lcg.h"
#include "openmc/ray.h"

using namespace openmc;

namespace {

// GeometryState sizes its coordinate levels from model::n_coord_levels, which
// is zero until a model is read. These tests do not trace through geometry,
// but they do construct rays, so one level has to exist.
struct CoordLevelFixture {
  CoordLevelFixture() : saved_(model::n_coord_levels)
  {
    model::n_coord_levels = 1;
  }
  ~CoordLevelFixture() { model::n_coord_levels = saved_; }
  int saved_;
};

ParticleRay make_ray(int64_t seed_id)
{
  return ParticleRay({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0},
    ParticleType::neutron(), 0.0, 1.0e6, seed_id);
}

} // namespace

TEST_CASE("ParticleRay initializes its RNG state")
{
  CoordLevelFixture fixture;

  ParticleRay ray = make_ray(0);

  // stream() indexes seeds_, so an uninitialized value is an out-of-bounds
  // pointer in current_seed(), which prn() writes through.
  REQUIRE(ray.stream() == STREAM_TRACKING);
  REQUIRE(ray.current_seed() == ray.seeds() + STREAM_TRACKING);

  // Every stream must be seeded, not just the tracking one: URR probability
  // tables draw from STREAM_URR_PTABLE during cross section lookup.
  uint64_t expected[N_STREAMS];
  init_particle_seeds(0, expected);
  for (int i = 0; i < N_STREAMS; ++i) {
    REQUIRE(ray.seeds(i) == expected[i]);
  }
}

TEST_CASE("ParticleRay seeding is deterministic and seed_id dependent")
{
  CoordLevelFixture fixture;

  ParticleRay a = make_ray(7);
  ParticleRay b = make_ray(7);
  ParticleRay c = make_ray(8);

  for (int i = 0; i < N_STREAMS; ++i) {
    REQUIRE(a.seeds(i) == b.seeds(i));
  }
  REQUIRE(a.seeds(STREAM_TRACKING) != c.seeds(STREAM_TRACKING));
}

TEST_CASE("ParticleRay copies its parent's RNG state without sharing it")
{
  CoordLevelFixture fixture;

  Particle parent;
  init_particle_seeds(1234, parent.seeds());
  parent.stream() = STREAM_TRACKING;
  parent.type() = ParticleType::neutron();
  parent.E() = 2.0e6;
  parent.time() = 3.0e-8;

  uint64_t parent_seeds_before[N_STREAMS];
  std::copy(parent.seeds(), parent.seeds() + N_STREAMS, parent_seeds_before);

  ParticleRay ray(parent, {1.0, 0.0, 0.0}, 2.0e6);

  // The ray starts from the parent's state, so the URR probability table
  // realization it samples matches the parent's.
  for (int i = 0; i < N_STREAMS; ++i) {
    REQUIRE(ray.seeds(i) == parent_seeds_before[i]);
  }
  REQUIRE(ray.stream() == STREAM_TRACKING);
  REQUIRE(ray.time() == parent.time());

  // Drawing from the ray must not disturb the parent's random walk.
  prn(ray.current_seed());
  for (int i = 0; i < N_STREAMS; ++i) {
    REQUIRE(parent.seeds(i) == parent_seeds_before[i]);
  }
  REQUIRE(ray.seeds(STREAM_TRACKING) != parent.seeds(STREAM_TRACKING));
}

TEST_CASE("ParticleRay accumulators start clean")
{
  CoordLevelFixture fixture;

  ParticleRay ray = make_ray(0);

  REQUIRE(ray.traversal_distance() == 0.0);
  REQUIRE(ray.traversal_mfp() == 0.0);
  REQUIRE_FALSE(ray.completed());

  // update_distance() in a void region must not touch the mfp.
  ray.update_distance(10.0);
  REQUIRE(ray.traversal_distance() == 10.0);
  REQUIRE(ray.traversal_mfp() == 0.0);
}
