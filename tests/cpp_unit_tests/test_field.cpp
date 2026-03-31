#include <string>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pugixml.hpp>

#include "openmc/mesh.h"
#include "openmc/field.h"

using namespace openmc;

TEST_CASE("Test TemperatureField functions with a regular mesh")
{
  // The XML data as a string
  std::string xml_string = R"(
        <mesh id="1">
            <dimension>2 2 2</dimension>
            <lower_left>-1 -1 -1</lower_left>
            <upper_right>1 1 1</upper_right>
      </mesh>
    )";

  // Create the mesh from a file
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_string(xml_string.c_str());
  pugi::xml_node root = doc.child("mesh");
  auto mesh = RegularMesh(root);

  // Define some temperature values
  vector<double> values = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0};

  // Create a temperature field
  TemperatureField temp_field = TemperatureField(&mesh, values);

  // Get temperature
  REQUIRE(temp_field.get_temperature(Position(0.5, 0.5, 0.5)) == 80.0);
  REQUIRE(temp_field.get_temperature(Position(-0.5, -0.5, -0.5)) == 10.0);
  REQUIRE(temp_field.get_temperature(Position(0.5, -0.5, -0.5)) == 20.0);
  REQUIRE(temp_field.get_temperature(Position(-0.5, -0.5, 0.5)) == 50.0);
  REQUIRE(temp_field.get_temperature(Position(0.0, 0.0, 0.0)) == 10.0);

  // Get sqrtkT
  REQUIRE(temp_field.get_sqrtkT(Position(0.5, 0.5, 0.5)) ==
          Catch::Approx(0.083029).margin(1.0E-6));

  // Update particle temperature
  Particle p;
  REQUIRE(p.sqrtkT() == -1.0);
  p.r() = Position(0.5, 0.5, 0.5);
  temp_field.update_particle_temperature(p);
  REQUIRE(p.sqrtkT() == Catch::Approx(0.083029).margin(1.0E-6));
  p.r() = Position(-0.5, -0.5, -0.5);
  temp_field.update_particle_temperature(p);
  REQUIRE(p.sqrtkT() == Catch::Approx(0.029355).margin(1.0E-6));
}
