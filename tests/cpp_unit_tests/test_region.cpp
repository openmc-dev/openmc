#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/particle_data.h"
#include "openmc/surface.h"

#include <fmt/core.h>
#include <pugixml.hpp>

#include <random>

namespace {

// Helper class to set up and tear down test surfaces
class SurfaceFixture {
public:
  SurfaceFixture()
  {
    pugi::xml_document doc;
    pugi::xml_node surf_node = doc.append_child("surface");
    surf_node.set_name("surface");
    surf_node.append_attribute("id") = "0";
    surf_node.append_attribute("type") = "x-plane";
    surf_node.append_attribute("coeffs") = "1";

    for (int i = 1; i < 10; ++i) {
      surf_node.attribute("id") = i;
      openmc::model::surfaces.push_back(
        std::make_unique<openmc::SurfaceXPlane>(surf_node));
      openmc::model::surface_map[i] = i - 1;
    }
  }

  ~SurfaceFixture()
  {
    openmc::model::surfaces.clear();
    openmc::model::surface_map.clear();
  }
};

// Helper class for testing multiple intersections with the same surface
class MultiIntersectionFixture {
public:
  MultiIntersectionFixture()
  {
    pugi::xml_document doc;

    auto plane = doc.append_child("surface");
    plane.append_attribute("id") = 1;
    plane.append_attribute("type") = "x-plane";
    plane.append_attribute("coeffs") = "5";
    openmc::model::surfaces.push_back(
      std::make_unique<openmc::SurfaceXPlane>(plane));
    openmc::model::surface_map[1] = 0;

    auto sphere = doc.append_child("surface");
    sphere.append_attribute("id") = 2;
    sphere.append_attribute("type") = "sphere";
    sphere.append_attribute("coeffs") = "5 0 0 1";
    openmc::model::surfaces.push_back(
      std::make_unique<openmc::SurfaceSphere>(sphere));
    openmc::model::surface_map[2] = 1;
  }

  ~MultiIntersectionFixture()
  {
    openmc::model::surfaces.clear();
    openmc::model::surface_map.clear();
  }
};

// Helper class for testing coincident surfaces with different representations
class CoincidentSurfaceFixture {
public:
  CoincidentSurfaceFixture()
  {
    pugi::xml_document doc;

    auto cylinder = doc.append_child("surface");
    cylinder.append_attribute("id") = 1;
    cylinder.append_attribute("type") = "z-cylinder";
    cylinder.append_attribute("coeffs") = "0 0 499.0000001";
    openmc::model::surfaces.push_back(
      std::make_unique<openmc::SurfaceZCylinder>(cylinder));
    openmc::model::surface_map[1] = 0;

    auto quadric = doc.append_child("surface");
    quadric.append_attribute("id") = 2;
    quadric.append_attribute("type") = "quadric";
    quadric.append_attribute("coeffs") =
      "1 1 7.498798913309288e-33 -7.498798913309288e-33 "
      "-1.2246467991473532e-16 -1.2246467991473532e-16 0 0 0 "
      "-249001.00009980003";
    openmc::model::surfaces.push_back(
      std::make_unique<openmc::SurfaceQuadric>(quadric));
    openmc::model::surface_map[2] = 1;
  }

  ~CoincidentSurfaceFixture()
  {
    openmc::model::surfaces.clear();
    openmc::model::surface_map.clear();
  }
};

// Helper class for testing a region outside of many finite rods
class RodsFixture {
public:
  static constexpr int N_RODS = 40;

  RodsFixture()
  {
    pugi::xml_document doc;
    auto add = [&](const char* type, std::string coeffs) {
      int id = openmc::model::surfaces.size() + 1;
      auto node = doc.append_child("surface");
      node.append_attribute("id") = id;
      node.append_attribute("type") = type;
      node.append_attribute("coeffs") = coeffs.c_str();
      if (type[0] == 'x') {
        openmc::model::surfaces.push_back(
          std::make_unique<openmc::SurfaceXPlane>(node));
      } else if (type[0] == 'y') {
        openmc::model::surfaces.push_back(
          std::make_unique<openmc::SurfaceYPlane>(node));
      } else if (type[2] == 'p') {
        openmc::model::surfaces.push_back(
          std::make_unique<openmc::SurfaceZPlane>(node));
      } else {
        openmc::model::surfaces.push_back(
          std::make_unique<openmc::SurfaceZCylinder>(node));
      }
      openmc::model::surface_map[id] = id - 1;
      return id;
    };

    // A box, and finite rods inside it, which may overlap
    region_spec = fmt::format("{} -{} {} -{} {} -{}", add("x-plane", "-10"),
      add("x-plane", "10"), add("y-plane", "-10"), add("y-plane", "10"),
      add("z-plane", "-10"), add("z-plane", "10"));
    std::mt19937 rng(1);
    std::uniform_real_distribution<double> center(-7.0, 7.0);
    for (int i = 0; i < N_RODS; ++i) {
      double x = center(rng);
      double y = center(rng);
      double z = center(rng);
      int cyl = add("z-cylinder", fmt::format("{} {} 1.5", x, y));
      int bottom = add("z-plane", fmt::format("{}", z - 2.0));
      int top = add("z-plane", fmt::format("{}", z + 2.0));
      region_spec += fmt::format(" ({} | -{} | {})", cyl, bottom, top);
    }
  }

  ~RodsFixture()
  {
    openmc::model::surfaces.clear();
    openmc::model::surface_map.clear();
  }

  std::string region_spec;
};

} // anonymous namespace

TEST_CASE("Test region simplification")
{
  SurfaceFixture fixture;

  SECTION("Original bug case from issue #3685")
  {
    // Input: "-1 2 (-3 4) | (-5 6)" was being incorrectly interpreted.
    // Nested intersections are merged and redundant parentheses are dropped.
    auto region = openmc::Region("(-1 2 (-3 4) | (-5 6))", 0);
    REQUIRE(region.str() == " ( -1 2 -3 4 ) | ( -5 6 )");
  }

  SECTION("Complement of a mixed expression")
  {
    // The complement applies to the grouped expression (1 2) | 3
    auto region = openmc::Region("~(1 2 | 3)", 0);
    REQUIRE(region.str() == " ( -1 | -2 ) -3");
  }

  SECTION("Complement of a parenthesized subexpression")
  {
    auto region = openmc::Region("4 ~(1 | 2 3)", 0);
    REQUIRE(region.str() == " 4 -1 ( -2 | -3 )");
  }

  SECTION("Simple union - no extra parentheses needed")
  {
    auto region = openmc::Region("1 | 2", 0);
    REQUIRE(region.str() == " 1 | 2");
  }

  SECTION("Intersection then union")
  {
    // Intersection should have higher precedence, so (1 2) grouped
    auto region = openmc::Region("1 2 | 3", 0);
    REQUIRE(region.str() == " ( 1 2 ) | 3");
  }

  SECTION("Union then intersection")
  {
    // The (2 3) intersection should be grouped
    auto region = openmc::Region("1 | 2 3", 0);
    REQUIRE(region.str() == " 1 | ( 2 3 )");
  }

  SECTION("Nested parentheses preserved")
  {
    // These parentheses are meaningful and should be preserved
    auto region = openmc::Region("(1 | 2) (3 | 4)", 0);
    REQUIRE(region.str() == " ( 1 | 2 ) ( 3 | 4 )");
  }

  SECTION("Deep nesting")
  {
    auto region = openmc::Region("((1 2) | (3 4)) 5", 0);
    REQUIRE(region.str() == " ( ( 1 2 ) | ( 3 4 ) ) 5");
  }

  SECTION("Multiple unions")
  {
    auto region = openmc::Region("1 | 2 | 3", 0);
    REQUIRE(region.str() == " 1 | 2 | 3");
  }

  SECTION("Multiple intersections")
  {
    auto region = openmc::Region("1 2 3", 0);
    // Simple cell - no operators in output
    REQUIRE(region.str() == " 1 2 3");
  }

  SECTION("Complex mixed expression")
  {
    auto region = openmc::Region("1 2 | 3 4 | 5 6", 0);
    REQUIRE(region.str() == " ( 1 2 ) | ( 3 4 ) | ( 5 6 )");
  }
}

TEST_CASE("Find boundary after virtual surface crossings")
{
  MultiIntersectionFixture fixture;
  openmc::Region region("-1 | -2", 0);

  SECTION("Starting inside the region")
  {
    // Along +x from x=1, entering the sphere at x=4 is virtual because x < 5.
    // Crossing the plane at x=5 is also virtual because the point is inside
    // the sphere. Exiting the sphere at x=6 is the first true boundary.
    auto [distance, surface] =
      region.distance({1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, 0);

    REQUIRE(distance == Catch::Approx(5.0));
    REQUIRE(surface == 2);
  }

  SECTION("Starting outside the region")
  {
    // Along -x from x=7, entering the sphere at x=6 is the first boundary.
    auto [distance, surface] =
      region.distance({7.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, 0);

    REQUIRE(distance == Catch::Approx(1.0));
    REQUIRE(surface == -2);
  }

  SECTION("Starting on a curved surface")
  {
    // Start on the sphere and travel obliquely through it. The plane crossing
    // is virtual, and accumulated roundoff must not cause the sphere exit to
    // be classified as another virtual crossing.
    auto [distance, surface] =
      region.distance({4.2, 0.6, 0.0}, {1.0, 0.0, 0.0}, -2);

    REQUIRE(distance == Catch::Approx(1.6));
    REQUIRE(surface == 2);
  }
}

TEST_CASE("Find boundary with and without working space")
{
  MultiIntersectionFixture fixture;
  openmc::Region region("-1 | -2", 0);

  // A particle normally has working space for every complex region. When it
  // does not, as can happen in event-based mode, the region is searched
  // without it and the boundary found is the same.
  openmc::GeometryState p;
  openmc::Position r[] = {{1.0, 0.0, 0.0}, {7.0, 0.0, 0.0}, {4.2, 0.6, 0.0}};
  openmc::Direction u[] = {{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
  int32_t on_surface[] = {0, 0, -2};
  for (int i = 0; i < 3; ++i) {
    auto expected = region.distance(r[i], u[i], on_surface[i]);
    p.surface_states().clear();
    auto without = region.distance(r[i], u[i], on_surface[i], &p);
    p.surface_states().resize(2);
    auto with = region.distance(r[i], u[i], on_surface[i], &p);
    REQUIRE(without.first == Catch::Approx(expected.first));
    REQUIRE(without.second == expected.second);
    REQUIRE(with.first == expected.first);
    REQUIRE(with.second == expected.second);
  }
}

TEST_CASE("Ignore roundoff-scale virtual surface crossings")
{
  CoincidentSurfaceFixture fixture;

  // These two surfaces describe effectively coincident cylinders, but the
  // small rotation in the general quadric causes their calculated
  // intersections to differ by roundoff. The nearby quadric intersection is
  // virtual and the next meaningful crossing is on the far side of the
  // cylinders.
  openmc::Region region("1 | -2", 0);
  auto [distance, surface] = region.distance(
    {-427.64056354508085, -257.1469395319449, -20.851278766740666},
    {0.8131471271523302, -0.39377148555937275, 0.4286441026822566}, -1);

  REQUIRE(distance == Catch::Approx(603.9161175466262));
  REQUIRE(surface == 1);
}

TEST_CASE("Find boundary of a region outside of many objects")
{
  // The region is an intersection of enough children with bounded boxes for
  // the boxes to be used. The boundary found along each ray must be where the
  // ray leaves the region.
  RodsFixture fixture;
  openmc::Region region(fixture.region_spec, 0);

  std::mt19937 rng(2);
  std::uniform_real_distribution<double> coord(-9.9, 9.9);
  std::uniform_real_distribution<double> mu(-1.0, 1.0);
  std::uniform_real_distribution<double> phi(0.0, 2.0 * openmc::PI);
  int n_inside = 0;
  for (int i = 0; i < 20000; ++i) {
    openmc::Position r {coord(rng), coord(rng), coord(rng)};
    double cos_theta = mu(rng);
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    double angle = phi(rng);
    openmc::Direction u {
      sin_theta * std::cos(angle), sin_theta * std::sin(angle), cos_theta};
    if (!region.contains(r, u, 0))
      continue;
    ++n_inside;

    // The region is inside a box, so a boundary is always found
    auto [distance, surface] = region.distance(r, u, 0);
    REQUIRE(distance < openmc::INFTY);
    REQUIRE(region.contains(r + (distance - 1e-6) * u, u, 0));
    REQUIRE(!region.contains(r + (distance + 1e-6) * u, u, 0));

    // The surface returned is the one crossed, with the sign of the side
    // entered
    openmc::Position hit = r + distance * u;
    const auto& surf = *openmc::model::surfaces[std::abs(surface) - 1];
    REQUIRE(std::abs(surf.evaluate(hit)) < 1e-6);
    REQUIRE((surface > 0) == (u.dot(surf.normal(hit)) > 0.0));
  }
  REQUIRE(n_inside > 1000);
}

TEST_CASE("Find boundary of a region outside of many objects with limited "
          "working space")
{
  // With working space for the surfaces of each child but not for all of the
  // surfaces, the children are searched separately when the ray starts in the
  // region, and the region is searched without working space otherwise. The
  // boundaries found are the same as with working space for all surfaces.
  RodsFixture fixture;
  openmc::Region region(fixture.region_spec, 0);
  openmc::GeometryState p;

  std::mt19937 rng(3);
  std::uniform_real_distribution<double> coord(-11.0, 11.0);
  std::uniform_real_distribution<double> mu(-1.0, 1.0);
  std::uniform_real_distribution<double> phi(0.0, 2.0 * openmc::PI);
  for (int i = 0; i < 2000; ++i) {
    openmc::Position r {coord(rng), coord(rng), coord(rng)};
    double cos_theta = mu(rng);
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    double angle = phi(rng);
    openmc::Direction u {
      sin_theta * std::cos(angle), sin_theta * std::sin(angle), cos_theta};

    auto expected = region.distance(r, u, 0);
    for (int n_states : {3, 0}) {
      p.surface_states().resize(n_states);
      auto [distance, surface] = region.distance(r, u, 0, &p);
      if (expected.first == openmc::INFTY) {
        REQUIRE(distance == openmc::INFTY);
      } else {
        REQUIRE(distance == Catch::Approx(expected.first));
        REQUIRE(surface == expected.second);
      }
    }
  }
}
