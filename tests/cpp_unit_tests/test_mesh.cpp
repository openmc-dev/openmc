#include <cstdio>
#include <iostream>
#include <map>
#include <string>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <pugixml.hpp>

#include "openmc/hdf5_interface.h"
#include "openmc/mesh.h"

#include "matchers/equals_position.hpp"

using namespace openmc;

TEST_CASE("Test mesh hdf5 roundtrip - regular")
{
  // The XML data as a string
  std::string xml_string = R"(
        <mesh id="1">
            <dimension>3 4 5</dimension>
            <lower_left>-2 -3 -5</lower_left>
            <upper_right>2 3 5</upper_right>
       </mesh>
    )";

  // Create a pugixml document object
  pugi::xml_document doc;

  // Load the XML from the string
  pugi::xml_parse_result result = doc.load_string(xml_string.c_str());

  pugi::xml_node root = doc.child("mesh");

  auto mesh = RegularMesh(root);

  hid_t file_id = file_open("mesh.h5", 'w');

  mesh.to_hdf5(file_id);

  file_close(file_id);

  hid_t file_id2 = file_open("mesh.h5", 'r');

  hid_t group = open_group(file_id2, "mesh 1");

  auto mesh2 = RegularMesh(group);

  file_close(file_id2);

  remove("mesh.h5");

  REQUIRE(mesh2.shape_ == mesh.shape_);

  REQUIRE(mesh2.lower_left() == mesh.lower_left());

  REQUIRE(mesh2.upper_right() == mesh.upper_right());
}

TEST_CASE("Test mesh hdf5 roundtrip - rectilinear")
{
  // The XML data as a string
  std::string xml_string = R"(
        <mesh id="1" type="rectilinear">
            <x_grid>0.0 1.0 5.0 10.0</x_grid>
            <y_grid>-10.0 -5.0 0.0</y_grid>
            <z_grid>-100.0 0.0 100.0</z_grid>
        </mesh>
    )";

  // Create a pugixml document object
  pugi::xml_document doc;

  // Load the XML from the string
  pugi::xml_parse_result result = doc.load_string(xml_string.c_str());

  pugi::xml_node root = doc.child("mesh");

  auto mesh = RectilinearMesh(root);

  hid_t file_id = file_open("mesh.h5", 'w');

  mesh.to_hdf5(file_id);

  file_close(file_id);

  hid_t file_id2 = file_open("mesh.h5", 'r');

  hid_t group = open_group(file_id2, "mesh 1");

  auto mesh2 = RectilinearMesh(group);

  file_close(file_id2);

  remove("mesh.h5");

  REQUIRE(mesh2.shape_ == mesh.shape_);

  REQUIRE(mesh2.grid_ == mesh.grid_);
}

TEST_CASE("Test mesh hdf5 roundtrip - cylindrical")
{
  // The XML data as a string
  std::string xml_string = R"(
        <mesh id="1" type="cylindrical">
            <r_grid>0.1 0.2 0.5 1.0</r_grid>
            <phi_grid>0.0 6.283185307179586</phi_grid>
            <z_grid>0.1 0.2 0.4 0.6 1.0</z_grid>
            <origin>0 0 0</origin>
        </mesh>
    )";

  // Create a pugixml document object
  pugi::xml_document doc;

  // Load the XML from the string
  pugi::xml_parse_result result = doc.load_string(xml_string.c_str());

  pugi::xml_node root = doc.child("mesh");

  auto mesh = CylindricalMesh(root);

  hid_t file_id = file_open("mesh.h5", 'w');

  mesh.to_hdf5(file_id);

  file_close(file_id);

  hid_t file_id2 = file_open("mesh.h5", 'r');

  hid_t group = open_group(file_id2, "mesh 1");

  auto mesh2 = CylindricalMesh(group);

  file_close(file_id2);

  remove("mesh.h5");

  REQUIRE(mesh2.shape_ == mesh.shape_);

  REQUIRE(mesh2.grid_ == mesh.grid_);
}

TEST_CASE("Test mesh hdf5 roundtrip - spherical")
{
  // The XML data as a string
  std::string xml_string = R"(
        <mesh id="1" type="spherical">
            <r_grid>0.1 0.2 0.5 1.0</r_grid>
            <theta_grid>0.0 3.141592653589793</theta_grid>
            <phi_grid>0.0 6.283185307179586</phi_grid>
            <origin>0.0 0.0 0.0</origin>
        </mesh>'
    )";

  // Create a pugixml document object
  pugi::xml_document doc;

  // Load the XML from the string
  pugi::xml_parse_result result = doc.load_string(xml_string.c_str());

  pugi::xml_node root = doc.child("mesh");

  auto mesh = SphericalMesh(root);

  hid_t file_id = file_open("mesh.h5", 'w');

  mesh.to_hdf5(file_id);

  file_close(file_id);

  hid_t file_id2 = file_open("mesh.h5", 'r');

  hid_t group = open_group(file_id2, "mesh 1");

  auto mesh2 = SphericalMesh(group);

  file_close(file_id2);

  remove("mesh.h5");

  REQUIRE(mesh2.shape_ == mesh.shape_);

  REQUIRE(mesh2.grid_ == mesh.grid_);
}

TEST_CASE("Test multiple meshes HDF5 roundtrip - spherical")
{
  // The XML data as a string
  std::string xml_string = R"(
        <meshes>
        <mesh id="1" type="spherical">
            <r_grid>0.1 0.2 0.5 1.0</r_grid>
            <theta_grid>0.0 3.141592653589793</theta_grid>
            <phi_grid>0.0 6.283185307179586</phi_grid>
            <origin>0.0 0.0 0.0</origin>
        </mesh>
         <mesh id="2">
            <dimension>3 4 5</dimension>
            <lower_left>-2 -3 -5</lower_left>
            <upper_right>2 3 5</upper_right>
       </mesh>
       </meshes>
    )";

  // Create a pugixml document object
  pugi::xml_document doc;

  // Load the XML from the string
  pugi::xml_parse_result result = doc.load_string(xml_string.c_str());

  pugi::xml_node root = doc.child("meshes");

  read_meshes(root);

  const auto spherical_mesh_xml =
    dynamic_cast<SphericalMesh*>(model::meshes[0].get());
  const auto regular_mesh_xml =
    dynamic_cast<RegularMesh*>(model::meshes[1].get());

  hid_t file_id = file_open("meshes.h5", 'w');

  hid_t root_group = create_group(file_id, "root");

  open_group(file_id, "root");

  meshes_to_hdf5(root_group);

  close_group(root_group);

  file_close(file_id);

  hid_t file_id2 = file_open("meshes.h5", 'r');

  hid_t root_group_read = open_group(file_id2, "root");

  hid_t mesh_group_read = open_group(root_group_read, "meshes");

  read_meshes(mesh_group_read);

  // increment mesh IDs to avoid collision during read
  for (auto& mesh : model::meshes) {
    mesh->set_id(mesh->id() + 10);
  }

  const auto spherical_mesh_hdf5 = dynamic_cast<SphericalMesh*>(
    model::meshes[model::mesh_map[spherical_mesh_xml->id_]].get());
  const auto regular_mesh_hdf5 = dynamic_cast<RegularMesh*>(
    model::meshes[model::mesh_map[regular_mesh_xml->id_]].get());

  remove("meshes.h5");

  REQUIRE(spherical_mesh_hdf5->shape_ == spherical_mesh_xml->shape_);
  REQUIRE(spherical_mesh_hdf5->grid_ == spherical_mesh_xml->grid_);

  REQUIRE(regular_mesh_hdf5->shape_ == regular_mesh_xml->shape_);
  REQUIRE(regular_mesh_hdf5->lower_left() == regular_mesh_xml->lower_left());
  REQUIRE(regular_mesh_hdf5->upper_right() == regular_mesh_xml->upper_right());
}

class RegularMeshFixture {
protected:
  RegularMesh mesh;

public:
  RegularMeshFixture()
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
    mesh = RegularMesh(root);

    // Add physical groups
    PGMap pg_map;
    vector<int> face_ids = {0, 24, 12, 36, 7, 31, 19, 43, 15, 39, 21, 45, 2, 26,
      8, 32, 4, 16, 10, 22, 29, 41, 35, 47};
    vector<int> physical_groups = {
      1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
    for (size_t i = 0; i < face_ids.size(); i++) {
      pg_map[physical_groups[i]].push_back(face_ids[i]);
    }
    mesh.pg_map() = pg_map;
  }
};

class RectilinearMeshFixture {
protected:
  RectilinearMesh mesh;

public:
  RectilinearMeshFixture()
  {
    // The XML data as a string
    std::string xml_string = R"(
          <mesh id="1" type="rectilinear">
              <x_grid>-1.0 0.0 1.0</x_grid>
              <y_grid>-1.0 0.0 1.0</y_grid>
              <z_grid>-1.0 0.0 1.0</z_grid>
          </mesh>
      )";

    // Create the mesh from a file
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_string(xml_string.c_str());
    pugi::xml_node root = doc.child("mesh");
    mesh = RectilinearMesh(root);
  }
};

#define GET_DISTANCE_TO_NEXT_BOUNDARY_TESTS                                    \
  int current_bin;                                                             \
  Position r;                                                                  \
  Position u;                                                                  \
  int next_bin;                                                                \
  double distance;                                                             \
                                                                               \
  SECTION("Test inside the mesh")                                              \
  {                                                                            \
    current_bin = 0;                                                           \
    r = Position(0.0, 0.0, 0.0);                                               \
    u = Position(1.0, 0.0, 0.0);                                               \
    distance = mesh.distance_to_next_boundary(current_bin, r, u, next_bin);    \
    REQUIRE(distance == 1.0);                                                  \
    REQUIRE(next_bin == -1);                                                   \
  }                                                                            \
                                                                               \
  SECTION("Test outside the mesh, going toward the mesh")                      \
  {                                                                            \
    current_bin = -1;                                                          \
    r = Position(-2.5, 0.0, 0.0);                                              \
    u = Position(1.0, 0.0, 0.0);                                               \
    distance = mesh.distance_to_next_boundary(current_bin, r, u, next_bin);    \
    REQUIRE(distance == 1.5);                                                  \
    REQUIRE(next_bin == 0);                                                    \
  }                                                                            \
                                                                               \
  SECTION("Test outside the mesh, not going toward the mesh")                  \
  {                                                                            \
    current_bin = -1;                                                          \
    r = Position(-2.0, 0.0, 0.0);                                              \
    u = Position(-1.0, 0.0, 0.0);                                              \
    distance = mesh.distance_to_next_boundary(current_bin, r, u, next_bin);    \
    REQUIRE(distance == INFTY);                                                \
    REQUIRE(next_bin == -1);                                                   \
  }                                                                            \
                                                                               \
  SECTION("Test on the mesh boundary, leaving the mesh")                       \
  {                                                                            \
    current_bin = 1;                                                           \
    r = Position(1.0, 0.0, 0.0);                                               \
    u = Position(1.0, 0.0, 0.0);                                               \
    distance = mesh.distance_to_next_boundary(current_bin, r, u, next_bin);    \
    REQUIRE(distance == INFTY);                                                \
    REQUIRE(next_bin == -1);                                                   \
  }                                                                            \
                                                                               \
  SECTION("Test close to the mesh boundary")                                   \
  {                                                                            \
    current_bin = 1;                                                           \
    r = Position(0.99999999999, 0.0, 0.0);                                     \
    u = Position(1.0, 0.0, 0.0);                                               \
    distance = mesh.distance_to_next_boundary(current_bin, r, u, next_bin);    \
    REQUIRE(distance == Catch::Approx(0.00000000001).margin(1.0E-12));         \
    REQUIRE(next_bin == -1);                                                   \
  }

TEST_CASE_METHOD(RegularMeshFixture, "Test distance_to_next_boundary()") {
  GET_DISTANCE_TO_NEXT_BOUNDARY_TESTS}

TEST_CASE_METHOD(RectilinearMeshFixture, "Test distance_to_next_boundary()") {
  GET_DISTANCE_TO_NEXT_BOUNDARY_TESTS}

#undef GET_DISTANCE_TO_NEXT_BOUNDARY_TESTS

TEST_CASE("Test get_index_in_direction - regular")
{
  // The XML data as a string
  std::string xml_string = R"(
        <mesh id="1">
            <dimension>2 1 1</dimension>
            <lower_left>-2 -2 -2</lower_left>
            <upper_right>2 2 2</upper_right>
       </mesh>
    )";

  // Create a pugixml document object
  pugi::xml_document doc;

  // Load the XML from the string
  pugi::xml_parse_result result = doc.load_string(xml_string.c_str());

  pugi::xml_node root = doc.child("mesh");

  auto mesh = RegularMesh(root);

  REQUIRE(mesh.get_index_in_direction(-4.1, 0) == 0);
  REQUIRE(mesh.get_index_in_direction(-4.0, 0) == 0);
  REQUIRE(mesh.get_index_in_direction(-3.9, 0) == 0);

  REQUIRE(mesh.get_index_in_direction(-2.1, 0) == 0);
  REQUIRE(mesh.get_index_in_direction(-2.0, 0) == 1); // lower left
  REQUIRE(mesh.get_index_in_direction(-1.9, 0) == 1);

  REQUIRE(mesh.get_index_in_direction(-0.1, 0) == 1);
  REQUIRE(mesh.get_index_in_direction(0.0, 0) == 1);
  REQUIRE(mesh.get_index_in_direction(0.1, 0) == 2);

  REQUIRE(mesh.get_index_in_direction(1.9, 0) == 2);
  REQUIRE(mesh.get_index_in_direction(2.0, 0) == 2); // upper right
  REQUIRE(mesh.get_index_in_direction(2.1, 0) == 3);

  REQUIRE(mesh.get_index_in_direction(3.9, 0) == 3);
  REQUIRE(mesh.get_index_in_direction(4.0, 0) == 3);
  REQUIRE(mesh.get_index_in_direction(4.1, 0) == 3);
}

TEST_CASE("Test get_index_in_direction - rectilinear")
{
  // The XML data as a string
  std::string xml_string = R"(
        <mesh id="1" type="rectilinear">
            <x_grid>-1.0 0.0 2.0</x_grid>
            <y_grid>-1.0 1.0</y_grid>
            <z_grid>-1.0 1.0</z_grid>
        </mesh>
    )";

  // Create a pugixml document object
  pugi::xml_document doc;

  // Load the XML from the string
  pugi::xml_parse_result result = doc.load_string(xml_string.c_str());

  pugi::xml_node root = doc.child("mesh");

  auto mesh = RectilinearMesh(root);

  REQUIRE(mesh.get_index_in_direction(-5.0, 0) == 0);

  REQUIRE(mesh.get_index_in_direction(-1.1, 0) == 0);
  REQUIRE(mesh.get_index_in_direction(-1.0, 0) == 1); // lower left
  REQUIRE(mesh.get_index_in_direction(-0.9, 0) == 1);

  REQUIRE(mesh.get_index_in_direction(1.9, 0) == 2);
  REQUIRE(mesh.get_index_in_direction(2.0, 0) == 2); // upper right
  REQUIRE(mesh.get_index_in_direction(2.1, 0) == 3);

  REQUIRE(mesh.get_index_in_direction(5.0, 0) == 3);
}

TEST_CASE("Test regular mesh ray tracing from outside")
{
  std::string xml_string = R"(
        <mesh id="1">
            <dimension>2 1 1</dimension>
            <lower_left>-2 -2 -2</lower_left>
            <upper_right>2 2 2</upper_right>
       </mesh>
    )";

  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_string(xml_string.c_str());
  auto mesh = RegularMesh(doc.child("mesh"));

  vector<int> bins;
  vector<double> lengths;

  mesh.bins_crossed(Position {-6.0, 0.0, 0.0}, Position {6.0, 0.0, 0.0},
    Direction {1.0, 0.0, 0.0}, bins, lengths);

  REQUIRE(bins == vector<int> {0, 1});
  REQUIRE(lengths.size() == 2);
  REQUIRE(lengths[0] == Catch::Approx(1.0 / 6.0));
  REQUIRE(lengths[1] == Catch::Approx(1.0 / 6.0));

  bins.clear();
  lengths.clear();

  mesh.bins_crossed(Position {6.0, 0.0, 0.0}, Position {-6.0, 0.0, 0.0},
    Direction {-1.0, 0.0, 0.0}, bins, lengths);

  REQUIRE(bins == vector<int> {1, 0});
  REQUIRE(lengths.size() == 2);
  REQUIRE(lengths[0] == Catch::Approx(1.0 / 6.0));
  REQUIRE(lengths[1] == Catch::Approx(1.0 / 6.0));
}

TEST_CASE_METHOD(RegularMeshFixture, "Test bins_and_surface_bins_crossed()")
{
  // Test cases:
  // - r0 inside and r1 inside
  // - r0 inside and r1 outside
  // - r0 outside and r1 inside
  // - r0 outside and r1 outside
  // - r0 inside and r1 inside: corner check
  auto [r0, r1, expected_leaving_surface_ids, expected_entering_surface_ids,
    expected_bins, expected_length_fractions] =
    GENERATE(table<Position, Position, vector<int>, vector<int>, vector<int>,
      vector<double>>({{Position(-0.5, -0.5, -0.5), Position(0.5, -0.5, -0.5),
                         {1}, {6}, {0, 1}, {0.5, 0.5}},
      {Position(-0.5, -0.5, -0.5), Position(2.0, -0.5, -0.5), {1, 7}, {6},
        {0, 1}, {0.2, 0.4}},
      {Position(-2.0, -0.5, -0.5), Position(0.5, -0.5, -0.5), {1}, {0, 6},
        {0, 1}, {0.4, 0.2}},
      {Position(-2.0, -0.5, -0.5), Position(2.0, -0.5, -0.5), {1, 7}, {0, 6},
        {0, 1}, {0.25, 0.25}},
      {Position(-0.5, -0.5, -0.5), Position(0.5, 0.5, 0.5), {1, 9, 23},
        {6, 20, 46}, {0, 1, 3, 7}, {0.5, 0.0, 0.0, 0.5}}}));

  vector<int> leaving_surface_ids;
  vector<int> entering_surface_ids;
  vector<int> bins;
  vector<double> length_fractions;

  mesh.bins_and_surface_bins_crossed(
    r0, r1, leaving_surface_ids, entering_surface_ids, bins, length_fractions);
  REQUIRE_THAT(
    leaving_surface_ids, Catch::Matchers::Equals(expected_leaving_surface_ids));
  REQUIRE_THAT(entering_surface_ids,
    Catch::Matchers::Equals(expected_entering_surface_ids));
  REQUIRE_THAT(bins, Catch::Matchers::Equals(expected_bins));
  REQUIRE(length_fractions.size() == expected_length_fractions.size());
  for (size_t i = 0; i < expected_length_fractions.size(); i++) {
    REQUIRE_THAT(length_fractions[i],
      Catch::Matchers::WithinAbs(expected_length_fractions[i], 1.0e-10));
  }
}

TEST_CASE_METHOD(RegularMeshFixture, "Test normalize_coordinates()")
{
  // r inside bin
  REQUIRE_THAT(mesh.normalize_coordinates(Position(-0.5, -0.5, -0.5), 0),
    EqualsPosition(Position(0.5, 0.5, 0.5), 1.0e-10));
  REQUIRE_THAT(mesh.normalize_coordinates(Position(-0.1, -0.7, -0.2), 0),
    EqualsPosition(Position(0.9, 0.3, 0.8), 1.0e-10));

  // r outside bin
  REQUIRE_THAT(mesh.normalize_coordinates(Position(0.5, -0.5, -0.5), 0),
    EqualsPosition(Position(1.0, 0.5, 0.5), 1.0e-10));
  REQUIRE_THAT(mesh.normalize_coordinates(Position(-0.5, -0.5, -0.5), 1),
    EqualsPosition(Position(0.0, 0.5, 0.5), 1.0e-10));
}

TEST_CASE_METHOD(RegularMeshFixture, "Test face_area()")
{
  REQUIRE(mesh.face_area(0) == 1.0);
  REQUIRE(mesh.face_area(1) == 1.0);
  REQUIRE(mesh.face_area(2) == 1.0);
  REQUIRE(mesh.face_area(3) == 1.0);
  REQUIRE(mesh.face_area(4) == 1.0);
  REQUIRE(mesh.face_area(5) == 1.0);

  // mesh.face_area(-1); // should fail
  // mesh.face_area(48); // should fail
}

TEST_CASE_METHOD(RegularMeshFixture, "Test sample_on_face()")
{
  // Check that sampled points are on the expected surface for different seeds

  uint64_t seed = GENERATE(1, 2);

  REQUIRE(mesh.sample_on_face(0, &seed).x == -1.0);
  REQUIRE(mesh.sample_on_face(1, &seed).x == 0.0);
  REQUIRE(mesh.sample_on_face(2, &seed).y == -1.0);
  REQUIRE(mesh.sample_on_face(3, &seed).y == 0.0);
  REQUIRE(mesh.sample_on_face(4, &seed).z == -1.0);
  REQUIRE(mesh.sample_on_face(5, &seed).z == 0.0);

  REQUIRE(mesh.sample_on_face(42, &seed).x == 0.0);
  REQUIRE(mesh.sample_on_face(43, &seed).x == 1.0);
  REQUIRE(mesh.sample_on_face(44, &seed).y == 0.0);
  REQUIRE(mesh.sample_on_face(45, &seed).y == 1.0);
  REQUIRE(mesh.sample_on_face(46, &seed).z == 0.0);
  REQUIRE(mesh.sample_on_face(47, &seed).z == 1.0);

  // mesh.sample_on_face(-1, &seed); // should fail
  // mesh.sample_on_face(100, &seed); // should fail
}

TEST_CASE_METHOD(RegularMeshFixture, "Regression test sample_on_face()")
{
  // Check that the other 2 coordinates are sampled

  uint64_t seed = 1;
  REQUIRE_THAT(mesh.sample_on_face(0, &seed),
    EqualsPosition(Position(-1.0, -0.2891826401, -0.2470039942), 1.0e-10));
}

TEST_CASE_METHOD(RegularMeshFixture, "Test return_vertex_unique_id()")
{
  // Lower left
  REQUIRE(mesh.return_vertex_unique_id({1, 1, 1}, 0) == 0);

  // Upper right
  REQUIRE(mesh.return_vertex_unique_id({2, 2, 2}, 7) == 26);

  // Middle from lower left bin
  REQUIRE(mesh.return_vertex_unique_id({1, 1, 1}, 7) == 13);

  // Middle from upper right bin
  REQUIRE(mesh.return_vertex_unique_id({2, 2, 2}, 0) == 13);

  // mesh.return_vertex_unique_id({0, 0, 0}, 0); // should fail
  // mesh.return_vertex_unique_id({3, 3, 3}, 0); // should fail
  // mesh.return_vertex_unique_id({1, 1, 1}, -1); // should fail
  // mesh.return_vertex_unique_id({1, 1, 1}, 8); // should fail
}

TEST_CASE_METHOD(RegularMeshFixture, "Test connectivity()")
{
  // Lower left element
  vector<int> expected = {0, 1, 3, 4, 9, 10, 12, 13};
  REQUIRE(mesh.connectivity(0) == expected);

  // Upper right element
  expected = {13, 14, 16, 17, 22, 23, 25, 26};
  REQUIRE(mesh.connectivity(7) == expected);

  // mesh.connectivity(-1); // should fail
  // mesh.connectivity(8); // should fail
}

#define GET_BIN_CLAMPED_TESTS                                                  \
  SECTION("r0 and r1 inside mesh - returns bin containing r1")                 \
  {                                                                            \
    Position r0 {-0.5, -0.5, -0.5};                                            \
    Position r1 {0.5, 0.5, -0.5};                                              \
    int bin0 = 0;                                                              \
    REQUIRE(mesh.get_bin_clamped(r0, r1, bin0) == 3);                          \
  }                                                                            \
                                                                               \
  SECTION(                                                                     \
    "r0 inside mesh and r1 outside mesh - returns last bin before leaving")    \
  {                                                                            \
    Position r0 {-0.5, -0.5, -0.5};                                            \
    Position r1 {1.5, 1.5, 1.5};                                               \
    int bin0 = 0;                                                              \
    REQUIRE(mesh.get_bin_clamped(r0, r1, bin0) == 7);                          \
  }                                                                            \
                                                                               \
  SECTION(                                                                     \
    "r0 inside mesh and r1 on internal boundary - returns bin containing r1")  \
  {                                                                            \
    Position r0 {-0.5, -0.5, -0.5};                                            \
    Position r1 {1.0, 1.0, 1.0};                                               \
    int bin0 = 0;                                                              \
    REQUIRE(mesh.get_bin_clamped(r0, r1, bin0) == 7);                          \
  }                                                                            \
                                                                               \
  SECTION("r1 == r0 - returns bin0 (no movement)")                             \
  {                                                                            \
    Position r0 {-0.5, -0.5, -0.5};                                            \
    Position r1 {-0.5, -0.5, -0.5};                                            \
    int bin0 = 0;                                                              \
    REQUIRE(mesh.get_bin_clamped(r0, r1, bin0) == bin0);                       \
  }                                                                            \
                                                                               \
  SECTION("r0 and r1 outside mesh with no traversal - returns bin0")           \
  {                                                                            \
    Position r0 {-1.5, -1.5, -1.5};                                            \
    Position r1 {-1.2, -1.2, -1.2};                                            \
    int bin0 = -1;                                                             \
    REQUIRE(mesh.get_bin_clamped(r0, r1, bin0) == bin0);                       \
  }                                                                            \
                                                                               \
  SECTION("r0 and r1 outside mesh with traversal - last bin before leaving")   \
  {                                                                            \
    Position r0 {-1.5, -1.5, -1.5};                                            \
    Position r1 {1.5, 1.5, 1.5};                                               \
    int bin0 = -1;                                                             \
    REQUIRE(mesh.get_bin_clamped(r0, r1, bin0) == 7);                          \
  }

TEST_CASE_METHOD(RegularMeshFixture, "Test get_bin_clamped()") {
  GET_BIN_CLAMPED_TESTS}

TEST_CASE_METHOD(RectilinearMeshFixture, "Test get_bin_clamped()") {
  GET_BIN_CLAMPED_TESTS}

#undef GET_BIN_CLAMPED_TESTS

// TODO: test get_bin_clamped with more mesh types

TEST_CASE_METHOD(RegularMeshFixture, "Test sample_on_physical_groups()")
{
  const double tol = 1.0e-12;
  uint64_t seed = 1;

  SECTION("Accept a single physical group")
  {
    // Using multiple calls of sample_on_physical_groups(), we verify that each
    // sample is on x=-1 (group 1)
    std::set<int> expected_bins = {0, 2, 4, 6};
    vector<int> physical_groups = {1};

    for (int i = 0; i < 1000; ++i) {
      Position p;
      int bin;
      mesh.sample_on_physical_groups(p, bin, &seed, physical_groups);

      REQUIRE_THAT(p.x, Catch::Matchers::WithinAbs(-1.0, tol));
      REQUIRE(expected_bins.count(bin) == 1);
      REQUIRE(p.y >= -1.0);
      REQUIRE(p.y <= 1.0);
      REQUIRE(p.z >= -1.0);
      REQUIRE(p.z <= 1.0);
    }
  }

  SECTION("Accept multiple physical groups")
  {
    // Using multiple calls of sample_on_physical_groups(), we verify that each
    // sample is either on x=-1 (group 1) or x=1 (group 2)
    vector<int> physical_groups = {1, 2};

    bool sampled_x_min = false;
    bool sampled_x_max = false;

    for (int i = 0; i < 1000; ++i) {
      Position p;
      int bin;
      mesh.sample_on_physical_groups(p, bin, &seed, physical_groups);

      bool on_x_min = std::abs(p.x - (-1.0)) < tol;
      bool on_x_max = std::abs(p.x - 1.0) < tol;
      REQUIRE((on_x_min || on_x_max));

      if (on_x_min) sampled_x_min = true;
      if (on_x_max) sampled_x_max = true;

      REQUIRE(p.y >= -1.0);
      REQUIRE(p.y <= 1.0);
      REQUIRE(p.z >= -1.0);
      REQUIRE(p.z <= 1.0);
    }

    // Both faces should have been sampled at least once
    REQUIRE(sampled_x_min);
    REQUIRE(sampled_x_max);
  }

  SECTION("Statistical test")
  {
    // Using multiple calls of sample_on_physical_groups(), we verify that each
    // of the 4 bins of physical group 1 gets ~25% of the samples, using a loose
    // tolerance (between 24% and 26%)
    vector<int> physical_groups = {1};

    std::map<int, int> bin_counts;
    const int n_samples = 10000;

    for (int i = 0; i < n_samples; ++i) {
      Position p;
      int bin;
      mesh.sample_on_physical_groups(p, bin, &seed, physical_groups);
      bin_counts[bin]++;
    }

    REQUIRE(bin_counts.size() == 4);
    for (auto& [bin, count] : bin_counts) {
      double fraction = static_cast<double>(count) / n_samples;
      REQUIRE(fraction > 0.24);
      REQUIRE(fraction < 0.26);
    }
  }
}
