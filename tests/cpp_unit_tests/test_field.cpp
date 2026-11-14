#include <string>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <pugixml.hpp>

#include "openmc/capi.h"
#include "openmc/field.h"
#include "openmc/geometry.h"
#include "openmc/mesh.h"
#include "openmc/nuclide.h"
#include "openmc/simulation.h"

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

  SECTION("Runtime temperature assignment through the C API")
  {
    // openmc_temperature_field_* act on the global field and consult the
    // global cross section temperature range, both of which have to be put
    // back however this section exits.
    struct GlobalStateGuard {
      double min {data::temperature_min};
      double max {data::temperature_max};
      ~GlobalStateGuard()
      {
        simulation::temperature_field = TemperatureField();
        data::temperature_min = min;
        data::temperature_max = max;
      }
    } guard;

    data::temperature_min = 300.0;
    data::temperature_max = 600.0;
    simulation::temperature_field = temp_field;

    double temperature;
    REQUIRE(openmc_temperature_field_set_temperature(0, 400.0) == 0);
    REQUIRE(openmc_temperature_field_get_temperature(0, &temperature) == 0);
    REQUIRE(temperature == 400.0);

    // Below the range covered by the loaded cross sections
    REQUIRE(openmc_temperature_field_set_temperature(0, 100.0) ==
            OPENMC_E_INVALID_ARGUMENT);

    // Above the range covered by the loaded cross sections
    REQUIRE(openmc_temperature_field_set_temperature(0, 1000.0) ==
            OPENMC_E_INVALID_ARGUMENT);

    // Not a finite, non-negative value
    REQUIRE(openmc_temperature_field_set_temperature(0, -1.0) ==
            OPENMC_E_INVALID_ARGUMENT);

    // Out of bounds index
    REQUIRE(openmc_temperature_field_set_temperature(8, 400.0) ==
            OPENMC_E_OUT_OF_BOUNDS);
  }

  SECTION("Distance to temperature mesh boundaries")
  {
    auto crossing = temp_field.next_mesh_crossing(
      0, Position {-0.5, -0.5, -0.5}, Direction {1.0, 0.0, 0.0});
    REQUIRE(crossing.distance == Catch::Approx(0.5));
    REQUIRE(crossing.next_bin == 1);

    crossing = temp_field.next_mesh_crossing(
      C_NONE, Position {-2.0, -0.5, -0.5}, Direction {1.0, 0.0, 0.0});
    REQUIRE(crossing.distance == Catch::Approx(1.0));
    REQUIRE(crossing.next_bin == 0);

    crossing = temp_field.next_mesh_crossing(
      0, Position {-1.0, -0.5, -0.5}, Direction {1.0, 0.0, 0.0});
    REQUIRE(crossing.distance == Catch::Approx(1.0));
    REQUIRE(crossing.next_bin == 1);

    // Reconcile stale state for a particle on the boundary entering the mesh.
    crossing = temp_field.next_mesh_crossing(
      C_NONE, Position {-1.0, -0.5, -0.5}, Direction {1.0, 0.0, 0.0});
    REQUIRE(crossing.distance == 0.0);
    REQUIRE(crossing.next_bin == 0);

    crossing = temp_field.next_mesh_crossing(
      1, Position {0.0, -0.5, -0.5}, Direction {1.0, 0.0, 0.0});
    REQUIRE(crossing.distance == Catch::Approx(1.0));
    REQUIRE(crossing.next_bin == C_NONE);

    constexpr double inverse_sqrt_three = 0.5773502691896258;
    crossing = temp_field.next_mesh_crossing(0, Position {-0.5, -0.5, -0.5},
      Direction {inverse_sqrt_three, inverse_sqrt_three, inverse_sqrt_three});
    REQUIRE(crossing.distance == Catch::Approx(0.5 / inverse_sqrt_three));
    REQUIRE(crossing.next_bin == 7);

    crossing = temp_field.next_mesh_crossing(
      C_NONE, Position {-2.0, 2.0, 0.0}, Direction {1.0, 0.0, 0.0});
    REQUIRE(crossing.distance == INFTY);
    REQUIRE(crossing.next_bin == C_NONE);
  }
}

TEST_CASE("Test temperature field crossings with a rectilinear mesh")
{
  std::string xml_string = R"(
        <mesh id="1" type="rectilinear">
            <x_grid>-1.0 0.0 2.0</x_grid>
            <y_grid>-1.0 1.0</y_grid>
            <z_grid>-1.0 1.0</z_grid>
        </mesh>
    )";

  pugi::xml_document doc;
  doc.load_string(xml_string.c_str());
  auto mesh = RectilinearMesh(doc.child("mesh"));
  TemperatureField field(&mesh, {300.0, 600.0});

  auto crossing = field.next_mesh_crossing(
    0, Position {-0.5, 0.0, 0.0}, Direction {1.0, 0.0, 0.0});
  REQUIRE(crossing.distance == Catch::Approx(0.5));
  REQUIRE(crossing.next_bin == 1);

  constexpr double inverse_sqrt_two = 0.7071067811865475;
  crossing = field.next_mesh_crossing(C_NONE, Position {-2.0, -2.0, 0.0},
    Direction {inverse_sqrt_two, inverse_sqrt_two, 0.0});
  REQUIRE(crossing.distance == Catch::Approx(1.0 / inverse_sqrt_two));
  REQUIRE(crossing.next_bin == 0);
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
    {// Negative temperature value -> error
      R"(
      <settings>
        <run_mode>eigenvalue</run_mode>
        <particles>200</particles>
        <batches>20</batches>
        <temperature_field>
          <mesh>1</mesh>
          <values>294.0 394.0 -1.0 594.0 694.0 794.0 894.0 994.0</values>
        </temperature_field>
        <mesh id="1">
          <dimension>2 2 2</dimension>
          <lower_left>0.0 0.0 0.0</lower_left>
          <upper_right>5.0 5.0 5.0</upper_right>
        </mesh>
      </settings>
      )",
      "Element 2 of the temperature field: Temperature of -1 K is not a "
      "finite, non-negative value."},
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
