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

TEST_CASE("Test Field - regular mesh", "[generators]")
{
  auto [mapping, values, outputs] =
    GENERATE(table<std::string, vector<double>, vector<double>>(
      {{"cell", {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0},
         {80.0, 10.0, 20.0, 50.0, 10.0, 10.0}},
        {"nodal",
          {10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
            21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0,
            32.0, 33.0, 34.0, 35.0, 36.0},
          {29.5, 16.5, 17.5, 25.5, 23.0, 10.0}}}));

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

  // Create a temperature field
  Field<double> field = Field<double>(&mesh, values, mapping);

  // Assign
  double saved_value = field.data().values()[2];
  field.assign(2, 150.0);
  REQUIRE(field.data().values()[2] == 150.0);
  field.data().values()[2] = saved_value; // Reassign back

  // Get bin
  REQUIRE(field.get_bin(Position(0.5, 0.5, 0.5)) == 7);
  REQUIRE(field.get_bin(Position(-0.5, -0.5, -0.5)) == 0);
  REQUIRE(field.get_bin(Position(0.5, -0.5, -0.5)) == 1);
  REQUIRE(field.get_bin(Position(-0.5, -0.5, 0.5)) == 4);
  REQUIRE(field.get_bin(Position(0.0, 0.0, 0.0)) == 0);
  REQUIRE(field.get_bin(Position(2.0, 2.0, 2.0)) == -1);

  // Evaluate (current position, current bin)
  REQUIRE(field.evaluate(Position(0.5, 0.5, 0.5), 7) == outputs[0]);
  REQUIRE(field.evaluate(Position(-0.5, -0.5, -0.5), 0) == outputs[1]);
  REQUIRE(field.evaluate(Position(0.5, -0.5, -0.5), 1) == outputs[2]);
  REQUIRE(field.evaluate(Position(-0.5, -0.5, 0.5), 4) == outputs[3]);
  REQUIRE(field.evaluate(Position(0.0, 0.0, 0.0), 0) == outputs[4]);

  // Trilinear interpolation
  if (mapping == "nodal") {
    REQUIRE(
      field.trilinear_interpolation(Position(0.5, 0.5, 0.5), 7) == outputs[0]);
    REQUIRE(field.trilinear_interpolation(Position(-0.5, -0.5, -0.5), 0) ==
            outputs[1]);
    REQUIRE(field.trilinear_interpolation(Position(0.5, -0.5, -0.5), 1) ==
            outputs[2]);
    REQUIRE(field.trilinear_interpolation(Position(-0.5, -0.5, 0.5), 4) ==
            outputs[3]);
    REQUIRE(
      field.trilinear_interpolation(Position(0.0, 0.0, 0.0), 0) == outputs[4]);
  }

  // Evaluate (previous position, current position, previous bin)
  REQUIRE(field.evaluate(Position(0.5, 0.5, 0.5), Position(-0.5, -0.5, -0.5),
            7) == outputs[1]);
  REQUIRE(field.evaluate(Position(-0.5, -0.5, -0.5), Position(0.5, 0.5, 0.5),
            0) == outputs[0]);
  REQUIRE(field.evaluate(Position(0.0, 0.0, 0.0), Position(-1.0, -1.0, -1.0),
            0) == outputs[5]);

  // Distance to next boundary
  int current_bin;
  Position r;
  Position u;
  int next_bin;
  double distance;

  // - Test inside the mesh
  current_bin = 0;
  r = Position(0.0, 0.0, 0.0);
  u = Position(1.0, 0.0, 0.0);
  distance = mesh.distance_to_next_boundary(current_bin, r, u, next_bin);
  REQUIRE(distance == 1.0);
  REQUIRE(next_bin == -1);

  // - Test outside the mesh, going toward the mesh
  current_bin = -1;
  r = Position(-2.5, 0.0, 0.0);
  u = Position(1.0, 0.0, 0.0);
  distance = mesh.distance_to_next_boundary(current_bin, r, u, next_bin);
  REQUIRE(distance == 1.5);
  REQUIRE(next_bin == 0);

  // - Test outside the mesh, not going toward the mesh
  current_bin = -1;
  r = Position(-2.0, 0.0, 0.0);
  u = Position(-1.0, 0.0, 0.0);
  distance = mesh.distance_to_next_boundary(current_bin, r, u, next_bin);
  REQUIRE(distance == INFTY);
  REQUIRE(next_bin == -1);

  // - Test on the mesh boundary, leaving the mesh
  current_bin = 1;
  r = Position(1.0, 0.0, 0.0);
  u = Position(1.0, 0.0, 0.0);
  distance = mesh.distance_to_next_boundary(current_bin, r, u, next_bin);
  REQUIRE(distance == INFTY);
  REQUIRE(next_bin == -1);

  // - Test close to the mesh boundary
  current_bin = 1;
  r = Position(0.99999999999, 0.0, 0.0);
  u = Position(1.0, 0.0, 0.0);
  distance = mesh.distance_to_next_boundary(current_bin, r, u, next_bin);
  REQUIRE(distance == Catch::Approx(0.00000000001).margin(1.0E-12));
  REQUIRE(next_bin == -1);
}

TEST_CASE("Test TemperatureField - regular mesh - cell-based only")
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
}

TEST_CASE("Test VelocityField - regular mesh")
{
  auto [mapping, values] = GENERATE(table<std::string, vector<Direction>>(
    {{"cell", {Direction(10.0, 10.0, 10.0), Direction(20.0, 20.0, 20.0),
                Direction(30.0, 30.0, 30.0), Direction(40.0, 40.0, 40.0),
                Direction(50.0, 50.0, 50.0), Direction(60.0, 60.0, 60.0),
                Direction(70.0, 70.0, 70.0), Direction(80.0, 80.0, 80.0)}},
      {"nodal", {Direction(10.0, 10.0, 10.0), Direction(11.0, 11.0, 11.0),
                  Direction(12.0, 12.0, 12.0), Direction(13.0, 13.0, 13.0),
                  Direction(14.0, 14.0, 14.0), Direction(15.0, 15.0, 15.0),
                  Direction(16.0, 16.0, 16.0), Direction(17.0, 17.0, 17.0),
                  Direction(18.0, 18.0, 18.0), Direction(19.0, 19.0, 19.0),
                  Direction(20.0, 20.0, 20.0), Direction(21.0, 21.0, 21.0),
                  Direction(22.0, 22.0, 22.0), Direction(23.0, 23.0, 23.0),
                  Direction(24.0, 24.0, 24.0), Direction(25.0, 25.0, 25.0),
                  Direction(26.0, 26.0, 26.0), Direction(27.0, 27.0, 27.0),
                  Direction(28.0, 28.0, 28.0), Direction(29.0, 29.0, 29.0),
                  Direction(30.0, 30.0, 30.0), Direction(31.0, 31.0, 31.0),
                  Direction(32.0, 32.0, 32.0), Direction(33.0, 33.0, 33.0),
                  Direction(34.0, 34.0, 34.0), Direction(35.0, 35.0, 35.0),
                  Direction(36.0, 36.0, 36.0)}}}));

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

  // Add physical group map
  mesh.pg_map() = {{1, {0, 24, 12, 36}}, {2, {7, 31, 19, 43}},
    {3, {15, 39, 21, 45, 2, 26, 8, 32, 4, 16}}, {4, {10, 22, 29, 41, 35, 47}}};

  // Create a temperature field
  VelocityField velocity_field = VelocityField(&mesh, values, mapping);

  // Add boundary conditions map
  velocity_field.bc_map() = {
    {BCType::INLET, {1}}, {BCType::OUTLET, {2}}, {BCType::WALL, {3, 4}}};

  // Get next bin
  Position r0;
  Position r1;
  int bin0;
  int bin1;
  BCType crossed_boundary;
  Position intersection;

  // - Next point is inside the mesh
  r0 = Position(-0.5, -0.5, -0.5);
  r1 = Position(0.5, 0.5, 0.5);
  bin0 = 0;
  bin1 =
    velocity_field.get_next_bin(r0, r1, bin0, crossed_boundary, intersection);
  REQUIRE(bin1 == 7);
  REQUIRE(crossed_boundary == BCType::NONE);
  REQUIRE(intersection == r1);

  // - Next point is outside the mesh
  r0 = Position(-0.5, -0.5, -0.5);
  r1 = Position(-0.5, -0.5, 1.5);
  bin0 = 0;
  bin1 =
    velocity_field.get_next_bin(r0, r1, bin0, crossed_boundary, intersection);
  REQUIRE(bin1 == -1);
  REQUIRE(crossed_boundary == BCType::WALL);
  REQUIRE(intersection == Position(-0.5, -0.5, 1.0));

  // Randomly place on inlet
  Position p = Position(0., 0., 0.);
  int bin = 0;
  uint64_t seed = 1;
  velocity_field.randomly_place_on_inlet(p, bin, &seed);
  REQUIRE(p.x == -1.0);
  REQUIRE(p.y == Catch::Approx(0.7529960058).margin(1.0E-10));
  REQUIRE(p.z == Catch::Approx(-0.6698810112).margin(1.0E-10));
  REQUIRE(bin == 2);

  // Get boundary conditions
  REQUIRE(velocity_field.get_boundary_condition(1) == BCType::INLET);
  REQUIRE(velocity_field.get_boundary_condition(2) == BCType::OUTLET);
  REQUIRE(velocity_field.get_boundary_condition(3) == BCType::WALL);
  REQUIRE(velocity_field.get_boundary_condition(4) == BCType::WALL);
}

TEST_CASE(
  "Test settings declaration exceptions - TemperatureField", "[generators]")
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
      "A mesh must be given for the temperature field."},
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
      "Temperature values must be given for the temperature field."},
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
