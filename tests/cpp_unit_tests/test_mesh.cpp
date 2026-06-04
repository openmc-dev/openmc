#include <cstdio>
#include <iostream>
#include <string>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <pugixml.hpp>

#include "openmc/hdf5_interface.h"
#include "openmc/mesh.h"

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
  RegularMeshFixture() {
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
  }
};

TEST_CASE_METHOD(RegularMeshFixture, "Test distance_to_next_boundary()")
{
  // Test distance_to_next_boundary()
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

class RectilinearMeshFixture {
protected:
  RectilinearMesh mesh;

public:
  RectilinearMeshFixture() {
    // The XML data as a string
    std::string xml_string = R"(
          <mesh id="1" type="rectilinear">
              <x_grid>-1.0 0.5 1.0</x_grid>
              <y_grid>-1.0 0.1 1.0</y_grid>
              <z_grid>-1.0 0.2 1.0</z_grid>
          </mesh>
      )";

    // Create the mesh from a file
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_string(xml_string.c_str());
    pugi::xml_node root = doc.child("mesh");
    mesh = RectilinearMesh(root);
  }
};

TEST_CASE_METHOD(RectilinearMeshFixture, "Test distance_to_next_boundary()")
{
  // Test distance_to_next_boundary()
  int current_bin;
  Position r;
  Position u;
  int next_bin;
  double distance;

  // - Test inside the mesh
  current_bin = 0;
  r = Position(0.5, 0.1, 0.2);
  u = Position(1.0, 0.0, 0.0);
  distance = mesh.distance_to_next_boundary(current_bin, r, u, next_bin);
  REQUIRE(distance == 0.5);
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

  Position result = mesh.sample_on_face(0, &seed);
  REQUIRE(result.x == -1.0);
  REQUIRE(result.y == Catch::Approx(-0.2891826401).margin(1.0E-10));
  REQUIRE(result.z == Catch::Approx(-0.2470039942).margin(1.0E-10));
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
