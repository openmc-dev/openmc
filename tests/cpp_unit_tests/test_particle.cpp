#include <catch2/catch_test_macros.hpp>

#include "openmc/geometry.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/particle_type.h"

namespace {

// Particle positions are stored in the coordinate levels, so give each particle
// one level while a test runs
class OneCoordLevel {
public:
  OneCoordLevel() : n_coord_levels_ {openmc::model::n_coord_levels}
  {
    openmc::model::n_coord_levels = 1;
  }
  ~OneCoordLevel() { openmc::model::n_coord_levels = n_coord_levels_; }

private:
  int n_coord_levels_;
};

void set_state(openmc::Particle& p, openmc::ParticleType type)
{
  p.type() = type;
  p.wgt() = 1.0;
  p.r() = {1.0, 2.0, 3.0};
  p.u() = {0.0, 0.0, 1.0};
  p.E() = 2.0e6;
  p.time() = 3.0e-6;
  p.lifetime() = 1.25e-6;
}

} // namespace

TEST_CASE("A source site starts with a zero lifetime")
{
  openmc::SourceSite site {};
  REQUIRE(site.lifetime == 0.0);
}

TEST_CASE_METHOD(
  OneCoordLevel, "Secondary and split neutrons keep the lifetime")
{
  using namespace openmc;

  Particle p;
  set_state(p, ParticleType::neutron());

  // A neutron from a reaction such as (n,2n)
  REQUIRE(p.create_secondary(p.wgt(), p.u(), p.E(), ParticleType::neutron()));
  // A neutron split by a weight window
  p.split(0.5);

  // Photons start their own clock
  REQUIRE(p.create_secondary(p.wgt(), p.u(), 1.0e6, ParticleType::photon()));
  p.type() = ParticleType::photon();
  p.split(0.5);

  // A neutron produced by a photon does not continue a neutron generation
  REQUIRE(p.create_secondary(p.wgt(), p.u(), 1.0e6, ParticleType::neutron()));

  const auto& bank = p.local_secondary_bank();
  REQUIRE(bank.size() == 5);
  for (int i = 0; i < 2; ++i) {
    REQUIRE(bank[i].particle == ParticleType::neutron());
    REQUIRE(bank[i].lifetime == 1.25e-6);
    REQUIRE(bank[i].time == 3.0e-6);
  }
  for (int i = 2; i < 4; ++i) {
    REQUIRE(bank[i].particle == ParticleType::photon());
    REQUIRE(bank[i].lifetime == 0.0);
    REQUIRE(bank[i].time == 3.0e-6);
  }
  REQUIRE(bank[4].particle == ParticleType::neutron());
  REQUIRE(bank[4].lifetime == 0.0);
  REQUIRE(bank[4].time == 3.0e-6);

  // A particle started from a banked site resumes the banked lifetime
  for (const auto& site : bank) {
    Particle q;
    q.lifetime() = -1.0;
    q.from_source(&site);
    REQUIRE(q.lifetime() == site.lifetime);
    REQUIRE(q.time() == site.time);
  }
}
