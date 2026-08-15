#include <string>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <pugixml.hpp>

#include "openmc/field.h"
#include "openmc/geometry.h"
#include "openmc/mesh.h"

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
  REQUIRE(temp_field.get_temperature(0) == 10.0);
  REQUIRE(temp_field.get_temperature(1) == 20.0);
  REQUIRE(temp_field.get_temperature(2) == 30.0);
  REQUIRE(temp_field.get_temperature(6) == 70.0);
  REQUIRE(temp_field.get_temperature(7) == 80.0);
  REQUIRE(temp_field.get_temperature(-1) == -1.0);

  // Get sqrtkT
  REQUIRE(temp_field.get_sqrtkT(7) == Catch::Approx(0.083029).margin(1.0E-6));

  // Get bin
  REQUIRE(temp_field.get_bin(Position(0.5, 0.5, 0.5)) == 7);
  REQUIRE(temp_field.get_bin(Position(-0.5, -0.5, -0.5)) == 0);
  REQUIRE(temp_field.get_bin(Position(0.5, -0.5, -0.5)) == 1);
  REQUIRE(temp_field.get_bin(Position(-0.5, -0.5, 0.5)) == 4);
  REQUIRE(temp_field.get_bin(Position(0.0, 0.0, 0.0)) == 0);
  REQUIRE(temp_field.get_bin(Position(2.0, 2.0, 2.0)) == -1);

  SECTION("Distance to temperature mesh boundaries")
  {
    int next_bin;

    auto distance = temp_field.distance_to_next_boundary(
      0, Position {-0.5, -0.5, -0.5}, Direction {1.0, 0.0, 0.0}, next_bin);
    REQUIRE(distance == Catch::Approx(0.5));
    REQUIRE(next_bin == 1);

    distance = temp_field.distance_to_next_boundary(
      -1, Position {-2.0, -0.5, -0.5}, Direction {1.0, 0.0, 0.0}, next_bin);
    REQUIRE(distance == Catch::Approx(1.0));
    REQUIRE(next_bin == 0);

    distance = temp_field.distance_to_next_boundary(
      -1, Position {-1.0, -0.5, -0.5}, Direction {1.0, 0.0, 0.0}, next_bin);
    REQUIRE(distance == 0.0);
    REQUIRE(next_bin == 0);

    distance = temp_field.distance_to_next_boundary(
      1, Position {0.0, -0.5, -0.5}, Direction {1.0, 0.0, 0.0}, next_bin);
    REQUIRE(distance == Catch::Approx(1.0));
    REQUIRE(next_bin == C_NONE);

    constexpr double inverse_sqrt_three = 0.5773502691896258;
    distance =
      temp_field.distance_to_next_boundary(0, Position {-0.5, -0.5, -0.5},
        Direction {inverse_sqrt_three, inverse_sqrt_three, inverse_sqrt_three},
        next_bin);
    REQUIRE(distance == Catch::Approx(0.5 / inverse_sqrt_three));
    REQUIRE(next_bin == 7);

    distance = temp_field.distance_to_next_boundary(
      -1, Position {-2.0, 2.0, 0.0}, Direction {1.0, 0.0, 0.0}, next_bin);
    REQUIRE(distance == INFTY);
    REQUIRE(next_bin == C_NONE);
  }
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
    {// Unsupported mesh type -> error
      R"(
      <settings>
        <run_mode>eigenvalue</run_mode>
        <particles>200</particles>
        <batches>20</batches>
        <temperature_field>
          <mesh>1</mesh>
          <values>294.0</values>
        </temperature_field>
        <mesh id="1" type="spherical">
          <r_grid>0.0 1.0</r_grid>
          <theta_grid>0.0 3.141592653589793</theta_grid>
          <phi_grid>0.0 6.283185307179586</phi_grid>
          <origin>0.0 0.0 0.0</origin>
        </mesh>
      </settings>
      )",
      "Temperature fields are only supported on regular and rectilinear "
      "meshes."},
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
