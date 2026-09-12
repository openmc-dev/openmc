#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "openmc/atomic_mass.h"
#include "openmc/constants.h"
#include "openmc/particle.h"
#include "openmc/particle_type.h"

using Catch::Approx;

TEST_CASE("Atomic and nuclear mass conventions")
{
  using namespace openmc;

  SECTION("light nuclides have AME2020 atomic masses")
  {
    REQUIRE(atomic_mass(1, 1) == Approx(1.007825031898));
    REQUIRE(atomic_mass(1, 2) == Approx(2.014101777844));
    REQUIRE(atomic_mass(1, 3) == Approx(3.01604928132));
    REQUIRE(atomic_mass(2, 3) == Approx(3.01602932197));
    REQUIRE(atomic_mass(2, 4) == Approx(4.00260325413));
  }

  SECTION("light nuclei have CODATA bare-particle masses")
  {
    REQUIRE(nuclear_mass(0, 1) == MASS_NEUTRON);
    REQUIRE(nuclear_mass(1, 1) == MASS_PROTON);
    REQUIRE(nuclear_mass(1, 2) == MASS_DEUTRON);
    REQUIRE(nuclear_mass(1, 3) == MASS_TRITON);
    REQUIRE(nuclear_mass(2, 3) == MASS_HELION);
    REQUIRE(nuclear_mass(2, 4) == MASS_ALPHA);
  }

  SECTION("heavier nuclear masses subtract the atomic electrons")
  {
    REQUIRE(
      nuclear_mass(26, 56) == Approx(atomic_mass(26, 56) - 26 * MASS_ELECTRON));
  }

  SECTION("invalid and unavailable masses are reported as missing")
  {
    REQUIRE(atomic_mass(0, 1) == 0.0);
    REQUIRE(atomic_mass(2, 1) == 0.0);
    REQUIRE(atomic_mass(60, 300) == 0.0);
    REQUIRE(nuclear_mass(2, 1) == 0.0);
    REQUIRE(nuclear_mass(60, 300) == 0.0);
  }
}

TEST_CASE("Particle masses are bare rest masses")
{
  using namespace openmc;

  REQUIRE(ParticleType::photon().mass() == 0.0);
  REQUIRE(ParticleType::electron().mass() == MASS_ELECTRON);
  REQUIRE(ParticleType::positron().mass() == MASS_ELECTRON);
  REQUIRE(ParticleType::neutron().mass() == MASS_NEUTRON);
  REQUIRE(ParticleType::proton().mass() == MASS_PROTON);
  REQUIRE(ParticleType::deuteron().mass() == MASS_DEUTRON);
  REQUIRE(ParticleType::triton().mass() == MASS_TRITON);
  REQUIRE(ParticleType::alpha().mass() == MASS_ALPHA);

  // Nuclear PDG codes refer to bare nuclei. The isomer digit is intentionally
  // ignored because ATOMIC_MASS contains ground-state masses only.
  REQUIRE(ParticleType {26, 56}.mass() == nuclear_mass(26, 56));
  REQUIRE(ParticleType {26, 56, 1}.mass() == nuclear_mass(26, 56));

  Particle p;
  p.type() = ParticleType {26, 56};
  REQUIRE(p.mass() == Approx(nuclear_mass(26, 56) * AMU_EV));
}
