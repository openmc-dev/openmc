#include <catch2/catch_test_macros.hpp>

#include "openmc/cell.h"
#include "openmc/surface.h"

#include <pugixml.hpp>

TEST_CASE("Test region simplification")
{
  pugi::xml_document doc;
  pugi::xml_node surf_node = doc.append_child("surface");
  surf_node.set_name("surface");
  surf_node.append_attribute("id") = "0";
  surf_node.append_attribute("type") = "x-plane";
  surf_node.append_attribute("coeffs") = "1";

  for (int i = 1; i < 7; ++i) {
    surf_node.attribute("id") = i;
    openmc::model::surfaces.push_back(
      std::make_unique<openmc::SurfaceXPlane>(surf_node));
    openmc::model::surface_map[i] = i - 1;
  }
  auto region = openmc::Region("(-1 2 (-3 4) | (-5 6))", 0);
  auto ref_val = " ( -1 2 -3 4 ) | ( -5 6 )";
  auto test_val = region.str();
  REQUIRE(test_val == ref_val);
}
