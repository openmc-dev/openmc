#include <string>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <pugixml.hpp>

#include "openmc/field.h"
#include "openmc/mesh.h"
#include "openmc/geometry.h"

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
  model::n_coord_levels = 1;
  Particle p;
  REQUIRE(p.sqrtkT() == -1.0);
  //p.r() = Position(0.5, 0.5, 0.5);
  //temp_field.update_particle_temperature(p);
  //REQUIRE(p.sqrtkT() == Catch::Approx(0.083029).margin(1.0E-6));
  //p.r() = Position(-0.5, -0.5, -0.5);
  //temp_field.update_particle_temperature(p);
  //REQUIRE(p.sqrtkT() == Catch::Approx(0.029355).margin(1.0E-6));
}

TEST_CASE("Test settings declaration exceptions for a temperature field",
  "[generators]")
{
  auto [input, error] = GENERATE(table<std::string, std::string>({
    {// If the number of values is not equal to the number of bins -> error
      R"(
      <settings>
        <run_mode>eigenvalue</run_mode>
        <particles>200</particles>
        <batches>20</batches>
        <temperature_field>
          <mesh>1</mesh>
          <values>294.0 394.0 494.0 594.0 694.0 794.0 894.0</values>
        </temperature_field>
        <mesh id="1">
          <dimension>2 2 2</dimension>
          <lower_left>0.0 0.0 0.0</lower_left>
          <upper_right>5.0 5.0 5.0</upper_right>
        </mesh>
      </settings>
      )",
      "Inconsistency in the temperature field: the number of values must be "
      "equal to the number of bins in the mesh."},
    {// If the mesh declared is not defined -> error
      R"(
      <settings>
        <run_mode>eigenvalue</run_mode>
        <particles>200</particles>
        <batches>20</batches>
        <temperature_field>
          <mesh>2</mesh>
          <values>294.0 394.0 494.0 594.0 694.0 794.0 894.0 994.0</values>
        </temperature_field>
        <mesh id="1">
          <dimension>2 2 2</dimension>
          <lower_left>0.0 0.0 0.0</lower_left>
          <upper_right>5.0 5.0 5.0</upper_right>
        </mesh>
      </settings>
      )",
      "Mesh 2 specified for the temperature field does not exist."},
    {// No mesh declared -> error
      R"(
      <settings>
        <run_mode>eigenvalue</run_mode>
        <particles>200</particles>
        <batches>20</batches>
        <temperature_field>
          <values>294.0 394.0 494.0 594.0 694.0 794.0 894.0 994.0</values>
        </temperature_field>
        <mesh id="1">
          <dimension>2 2 2</dimension>
          <lower_left>0.0 0.0 0.0</lower_left>
          <upper_right>5.0 5.0 5.0</upper_right>
        </mesh>
      </settings>
      )",
      "A mesh should be given for the temperature field."},
    {// No values declared -> error
      R"(
      <settings>
        <run_mode>eigenvalue</run_mode>
        <particles>200</particles>
        <batches>20</batches>
        <temperature_field>
          <mesh>1</mesh>
        </temperature_field>
        <mesh id="1">
          <dimension>2 2 2</dimension>
          <lower_left>0.0 0.0 0.0</lower_left>
          <upper_right>5.0 5.0 5.0</upper_right>
        </mesh>
      </settings>
      )",
      "Temperature values should be given for the temperature field."},
  }));

  free_memory_mesh();
  free_memory_settings();
  settings::run_mode = RunMode::UNSET;

  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_string(input.c_str());
  pugi::xml_node root = doc.child("settings");

  CAPTURE(input);
  REQUIRE_THROWS_WITH(read_settings_xml(root), error);
  doc.reset();
}
