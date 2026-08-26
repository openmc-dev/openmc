#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "openmc/dnp_drift.h"
#include "openmc/mesh.h"
#include "openmc/position.h"
#include "openmc/simulation.h"
#include "openmc/streamline_integrator.h"

using namespace openmc;

// ----------------------------------------------------------------------------
// _adjust_position
// ----------------------------------------------------------------------------

TEST_CASE("Test linear position adjustment")
{
  // Typical use case
  Position a = Position(-1.0, -1.0, -1.0);
  Position b = Position(1.0, 1.0, 1.0);
  _adjust_position(a, b, 3.0, 1.0, 2.5);
  REQUIRE(a == Position(0.0, 0.0, 0.0));
}

// ----------------------------------------------------------------------------
// _adjust_time
// ----------------------------------------------------------------------------

TEST_CASE("Test linear time adjustment")
{
  // Typical use case
  double t = 3.0;
  Position a = Position(-1.0, -1.0, -1.0);
  Position b = Position(1.0, 1.0, 1.0);
  Position c = Position(0.0, 0.0, 0.0);
  _adjust_time(t, 2.0, a, b, c);
  REQUIRE(t == 2.5);
}

// ----------------------------------------------------------------------------
// transport_dnp
// ----------------------------------------------------------------------------

class DNPTransportModelFixture {
protected:
  RegularMesh mesh;
  BCMap bc_map;
  SourceSite site;
  uint64_t seed = 1;

  void set_site(const Position& r, const Direction& u)
  {
    site.r = r;
    site.u = u;
  }

  void set_velocity_field(const Direction& v, const std::string& type = "cell")
  {
    vector<Direction> velocities(8, v);
    simulation::velocity_field = VelocityField(&mesh, velocities, type);
    simulation::velocity_field.bc_map() = bc_map;
  }

  void set_integrator(double step_size)
  {
    simulation::streamline_integrator =
      std::make_unique<RK4StreamlineIntegrator>(step_size);
  }

public:
  DNPTransportModelFixture()
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

    // Boundary condition map
    bc_map[BCType::INLET] = {1};
    bc_map[BCType::OUTLET] = {2};
    bc_map[BCType::WALL] = {3, 4};

    // Default simulation setup
    simulation::streamline_integrator =
      std::make_unique<RK4StreamlineIntegrator>(0.8);
    set_velocity_field(Direction(1.0, 0.0, 0.0));

    // Default site
    set_site(Position(0.0, 0.0, 0.0), Direction(1.0, 0.0, 0.0));
  }

  ~DNPTransportModelFixture()
  {
    // Reset all openmc global variables
    simulation::streamline_integrator.reset();
    simulation::velocity_field = VelocityField();
    settings::dnp_drift_external_travel_time = 0.0;
    settings::dnp_drift_recycling_on = false;
  }
};

TEST_CASE_METHOD(DNPTransportModelFixture,
  "Test DNP transport - cross outlet with recycling "
  "on - decay inside the system")
{
  settings::dnp_drift_external_travel_time = 1.0;
  settings::dnp_drift_recycling_on = true;

  bool is_inside = transport_dnp(site, 2.5, &seed);

  REQUIRE(is_inside);
  REQUIRE(site.r.x == Catch::Approx(-0.5).margin(1.0E-10));
}

TEST_CASE_METHOD(DNPTransportModelFixture,
  "Test DNP transport - cross outlet with recycling "
  "on - decay outside the system")
{
  settings::dnp_drift_external_travel_time = 100.0;
  settings::dnp_drift_recycling_on = true;

  bool is_inside = transport_dnp(site, 2.5, &seed);

  REQUIRE(!is_inside);
  REQUIRE(site.r == Position(1.0, 0.0, 0.0));
}

TEST_CASE_METHOD(DNPTransportModelFixture,
  "Test DNP transport - cross outlet with recycling off")
{
  bool is_inside = transport_dnp(site, 2.5, &seed);

  REQUIRE(!is_inside);
  REQUIRE(site.r == Position(1.0, 0.0, 0.0));
}

TEST_CASE_METHOD(DNPTransportModelFixture, "Test DNP transport - cross wall")
{
  set_velocity_field(Direction(0.0, 1.0, 0.0));

  bool is_inside = transport_dnp(site, 10.0, &seed);

  REQUIRE(is_inside);
  REQUIRE(site.r == Position(0.0, 1.0, 0.0));
}

TEST_CASE_METHOD(DNPTransportModelFixture, "Test DNP transport - cross inlet")
{
  simulation::velocity_field.assign(0, Direction(-1.0, 0.0, 0.0));
  set_site(Position(-0.5, -0.5, -0.5), Direction(1.0, 0.0, 0.0));

  bool is_inside = transport_dnp(site, 1.0, &seed);

  REQUIRE(is_inside);
  REQUIRE(site.r.x == Catch::Approx(-0.5).margin(1.0E-10));
}

TEST_CASE_METHOD(DNPTransportModelFixture,
  "Test DNP transport - stop exactly on mesh boundary")
{
  settings::dnp_drift_external_travel_time = 1.0;
  settings::dnp_drift_recycling_on = true;
  set_integrator(1.0);

  bool is_inside = transport_dnp(site, 3.0, &seed);

  REQUIRE(is_inside);
  REQUIRE(site.r.x == Catch::Approx(0.0).margin(1.0E-10));
}

// ----------------------------------------------------------------------------
// reconcile_precursor_drift
// ----------------------------------------------------------------------------

// TODO: prepare model for exhaustive_find_cell() to work
