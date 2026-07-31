#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "openmc/cell.h"
#include "openmc/surface.h"

#include <pugixml.hpp>

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

} // anonymous namespace

TEST_CASE("Test region simplification")
{
  SurfaceFixture fixture;

  SECTION("Original bug case from issue #3685")
  {
    // Input: "-1 2 (-3 4) | (-5 6)" was being incorrectly interpreted
    auto region = openmc::Region("(-1 2 (-3 4) | (-5 6))", 0);
    REQUIRE(region.str() == " ( ( -1 2 ( -3 4 ) ) | ( -5 6 ) )");
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
