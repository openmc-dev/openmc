#include <cmath>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <pugixml.hpp>

#include "openmc/bank.h"
#include "openmc/cell.h"
#include "openmc/event.h"
#include "openmc/field.h"
#include "openmc/geometry.h"
#include "openmc/geometry_aux.h"
#include "openmc/mesh.h"
#include "openmc/particle.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/surface.h"

using namespace openmc;

namespace {

class CountingMesh : public RegularMesh {
public:
  using RegularMesh::RegularMesh;
  int get_bin(Position r) const override
  {
    ++lookups;
    return RegularMesh::get_bin(r);
  }
  mutable int lookups {0};
};

class FieldTransportFixture {
public:
  FieldTransportFixture(const char* geometry_xml = R"(<geometry>
      <surface id="1" type="x-plane" coeffs="-2" boundary="vacuum"/>
      <surface id="2" type="x-plane" coeffs="2" boundary="vacuum"/>
      <cell id="1" material="void" region="1 -2" temperature="250"/>
    </geometry>)")
  {
    settings::run_mode = RunMode::FIXED_SOURCE;
    settings::weight_window_checkpoint_surface = false;
    settings::use_shared_secondary_bank = false;
    settings::time_cutoff.fill(INFTY);
    pugi::xml_document geometry;
    geometry.load_string(geometry_xml);
    read_geometry_xml(geometry.document_element());
    model::universe_level_counts.clear();
    finalize_geometry();
    finalize_cell_densities();

    pugi::xml_document doc;
    doc.load_string(R"(<mesh id="1">
      <dimension>2 1 1</dimension><lower_left>0 -1 -1</lower_left>
      <upper_right>2 1 1</upper_right></mesh>)");
    mesh = make_unique<CountingMesh>(doc.document_element());
    simulation::temperature_field = TemperatureField(mesh.get(), {300, 600});
    settings::temperature_field_on = true;
  }

  ~FieldTransportFixture()
  {
    free_event_queues();
    simulation::temperature_field = TemperatureField();
    free_memory_geometry();
    free_memory_surfaces();
    model::universe_level_counts.clear();
    settings::temperature_field_on = field_on_;
    settings::run_mode = run_mode_;
    settings::time_cutoff = cutoff_;
    settings::weight_window_checkpoint_surface = ww_;
    settings::use_shared_secondary_bank = shared_;
  }

  void initialize(
    Particle& p, Position r = {0.25, 0, 0}, Direction u = {1, 0, 0})
  {
    p.init_from_r_u(r, u);
    p.type() = ParticleType::neutron();
    p.E() = 1.e6;
    p.wgt() = 1;
    p.time() = 0;
    p.stream() = STREAM_TRACKING;
    p.seeds(0) = 123;
    REQUIRE(exhaustive_find_cell(p));
    p.macro_xs().total = 0;
    p.macro_xs().nu_fission = 0;
  }

  void collision_at(Particle& p, double distance)
  {
    auto seed = *p.current_seed();
    p.macro_xs().total = -std::log(prn(&seed)) / distance;
  }

  unique_ptr<CountingMesh> mesh;

private:
  bool field_on_ {settings::temperature_field_on};
  RunMode run_mode_ {settings::run_mode};
  decltype(settings::time_cutoff) cutoff_ {settings::time_cutoff};
  bool ww_ {settings::weight_window_checkpoint_surface};
  bool shared_ {settings::use_shared_secondary_bank};
};

} // namespace

TEST_CASE("Select the actual event after comparing transport distances")
{
  FieldTransportFixture fixture;
  Particle p;
  fixture.initialize(p);

  SECTION("Mesh crossing")
  {
    auto result = p.event_advance();
    REQUIRE(result.event == TransportEvent::TEMPERATURE_FIELD);
    REQUIRE(result.next_temperature_field_bin == 1);
    REQUIRE(p.r().x == Catch::Approx(1));
  }
  SECTION("Geometry without a field")
  {
    settings::temperature_field_on = false;
    REQUIRE(p.event_advance().event == TransportEvent::GEOMETRY);
    REQUIRE(p.r().x == Catch::Approx(2));
  }
  SECTION("Collision before the candidate mesh crossing")
  {
    fixture.collision_at(p, 0.25);
    REQUIRE(p.event_advance().event == TransportEvent::COLLISION);
    REQUIRE(p.r().x == Catch::Approx(0.5));
  }
  SECTION("Collision wins a tie with a mesh crossing")
  {
    fixture.collision_at(p, 0.75);
    REQUIRE(p.event_advance().event == TransportEvent::COLLISION);
  }
  SECTION("Cutoff wins a tie with a mesh crossing")
  {
    settings::time_cutoff[0] = 0.75 / p.speed();
    REQUIRE(p.event_advance().event == TransportEvent::TIME_CUTOFF);
    REQUIRE_FALSE(p.alive());
  }
}

TEST_CASE("An event-mode cutoff revives a local secondary")
{
  FieldTransportFixture fixture;
  init_event_queues(1);
  auto& p = simulation::particles[0];
  fixture.initialize(p);
  settings::time_cutoff[0] = 0.25 / p.speed();
  SourceSite secondary;
  secondary.r = {0.5, 0, 0};
  secondary.u = {1, 0, 0};
  secondary.E = 1.e6;
  secondary.wgt = 1;
  secondary.particle = ParticleType::neutron();
  p.local_secondary_bank().push_back(secondary);
  simulation::advance_particle_queue.thread_safe_append({p, 0});
  process_advance_particle_events();
  REQUIRE(p.alive());
  REQUIRE(p.local_secondary_bank().empty());
  REQUIRE(simulation::calculate_nonfuel_xs_queue.size() == 1);
  REQUIRE(simulation::collision_queue.size() == 0);
  REQUIRE(simulation::surface_crossing_queue.size() == 0);
}

TEST_CASE(
  "Geometry transmission reuses bins and explicitly crosses coincident meshes")
{
  const double split = GENERATE(-1.0, 0.0, 0.5, 1.0, 2.0);
  const bool event_based = GENERATE(false, true);
  CAPTURE(split, event_based);
  auto xml = fmt::format(R"(<geometry>
    <surface id="1" type="x-plane" coeffs="-2" boundary="vacuum"/>
    <surface id="2" type="x-plane" coeffs="3" boundary="vacuum"/>
    <surface id="3" type="x-plane" coeffs="{}"/>
    <cell id="1" material="void" region="1 -3" temperature="250"/>
    <cell id="2" material="void" region="3 -2" temperature="275"/>
  </geometry>)",
    split);
  FieldTransportFixture fixture(xml.c_str());
  init_event_queues(1);
  auto& p = simulation::particles[0];
  fixture.initialize(p, {split - 0.25, 0, 0});
  const double previous = p.sqrtkT();
  fixture.mesh->lookups = 0;
  const int expected_bin = split < 0 || split == 2 ? C_NONE : split < 1 ? 0 : 1;
  const double temperature = expected_bin == C_NONE ? 275
                             : expected_bin == 0    ? 300
                                                    : 600;

  if (event_based) {
    simulation::advance_particle_queue.thread_safe_append({p, 0});
    process_advance_particle_events();
    REQUIRE(simulation::surface_crossing_queue.size() == 1);
    process_surface_crossing_events();
    REQUIRE(p.n_event() == 1);
    REQUIRE(simulation::calculate_nonfuel_xs_queue.size() == 1);
  } else {
    auto result = p.event_advance();
    REQUIRE(result.event == TransportEvent::GEOMETRY);
    REQUIRE(result.temperature_field_crossing ==
            (split == 0 || split == 1 || split == 2));
    p.event_cross_surface(result);
  }
  REQUIRE(p.lowest_coord().cell() == 1);
  REQUIRE(p.temperature_field_bin() == expected_bin);
  REQUIRE(p.sqrtkT() == Catch::Approx(std::sqrt(K_BOLTZMANN * temperature)));
  REQUIRE(p.sqrtkT_last() == previous);
  REQUIRE(fixture.mesh->lookups == 0);
}

TEST_CASE("Coincident boundary conditions resolve the outgoing field bin")
{
  const auto bc = GENERATE("reflective", "white", "periodic", "vacuum");
  CAPTURE(bc);
  auto xml = fmt::format(R"(<geometry>
    <surface id="1" type="x-plane" coeffs="0" boundary="{}"/>
    <surface id="2" type="x-plane" coeffs="2" boundary="{}"/>
    <cell id="1" material="void" region="1 -2" temperature="250"/>
  </geometry>)",
    bc, bc);
  FieldTransportFixture fixture(xml.c_str());
  Particle p;
  fixture.initialize(p, {1.25, 0, 0});
  fixture.mesh->lookups = 0;
  auto result = p.event_advance();
  REQUIRE(result.event == TransportEvent::GEOMETRY);
  REQUIRE(result.temperature_field_crossing);
  REQUIRE(result.next_temperature_field_bin == C_NONE);
  p.event_cross_surface(result);
  if (std::string(bc) == "vacuum") {
    REQUIRE_FALSE(p.alive());
    REQUIRE(fixture.mesh->lookups == 0);
  } else {
    const int bin = std::string(bc) == "periodic" ? 0 : 1;
    REQUIRE(p.temperature_field_bin() == bin);
    REQUIRE(p.sqrtkT() == simulation::temperature_field.get_sqrtkT(bin));
    REQUIRE(fixture.mesh->lookups == 1);
  }
}

TEST_CASE("Reflection at a mesh corner uses the outgoing direction")
{
  FieldTransportFixture fixture(R"(<geometry>
    <surface id="1" type="x-plane" coeffs="-2" boundary="vacuum"/>
    <surface id="2" type="x-plane" coeffs="1" boundary="reflective"/>
    <cell id="1" material="void" region="1 -2" temperature="250"/>
  </geometry>)");
  pugi::xml_document doc;
  doc.load_string(R"(<mesh id="2"><dimension>2 2 1</dimension>
    <lower_left>0 -1 -1</lower_left><upper_right>2 1 1</upper_right></mesh>)");
  fixture.mesh = make_unique<CountingMesh>(doc.document_element());
  simulation::temperature_field =
    TemperatureField(fixture.mesh.get(), {300, 400, 500, 600});
  Particle p;
  const double u = std::sqrt(0.5);
  fixture.initialize(p, {0.5, -0.5, 0}, {u, u, 0});
  auto result = p.event_advance();
  REQUIRE(result.next_temperature_field_bin == 3);
  p.event_cross_surface(result);
  REQUIRE(p.temperature_field_bin() == 2);
  REQUIRE(p.sqrtkT() == simulation::temperature_field.get_sqrtkT(2));
}

TEST_CASE("General geometry queries locate the temperature field")
{
  FieldTransportFixture fixture;
  Particle p;
  fixture.initialize(p);
  p.r().x = 1.25;
  fixture.mesh->lookups = 0;
  REQUIRE(exhaustive_find_cell(p));
  REQUIRE(p.temperature_field_bin() == 1);
  REQUIRE(fixture.mesh->lookups == 1);
  p.r().x = -0.25;
  REQUIRE(exhaustive_find_cell(p));
  REQUIRE(p.temperature_field_bin() == C_NONE);
  REQUIRE(p.sqrtkT() == Catch::Approx(std::sqrt(K_BOLTZMANN * 250)));
}

TEST_CASE("A stale outside status is reconciled once without a geometry lookup")
{
  FieldTransportFixture fixture;
  Particle p;
  fixture.initialize(p, {0, 0, 0});
  p.temperature_field_bin() = C_NONE;
  fixture.mesh->lookups = 0;
  auto result = p.event_advance();
  REQUIRE(result.event == TransportEvent::TEMPERATURE_FIELD);
  REQUIRE(p.boundary().distance() == 0);
  p.event_cross_temperature_field(result.next_temperature_field_bin);
  REQUIRE(p.temperature_field_bin() == 0);
  REQUIRE(p.event_advance().event == TransportEvent::TEMPERATURE_FIELD);
  REQUIRE(p.boundary().distance() == Catch::Approx(1));
  REQUIRE(fixture.mesh->lookups == 0);
}

TEST_CASE(
  "A periodic mapping relocates even when it precedes the mesh boundary")
{
  FieldTransportFixture fixture(R"(<geometry>
    <surface id="1" type="x-plane" coeffs="0" boundary="periodic"/>
    <surface id="2" type="x-plane" coeffs="2" boundary="periodic"/>
    <cell id="1" material="void" region="1 -2" temperature="250"/>
  </geometry>)");
  pugi::xml_document doc;
  doc.load_string(R"(<mesh id="2"><dimension>2 1 1</dimension>
    <lower_left>-1 -1 -1</lower_left><upper_right>3 1 1</upper_right></mesh>)");
  fixture.mesh = make_unique<CountingMesh>(doc.document_element());
  simulation::temperature_field =
    TemperatureField(fixture.mesh.get(), {300, 600});
  Particle p;
  fixture.initialize(p, {1.25, 0, 0});
  fixture.mesh->lookups = 0;
  auto result = p.event_advance();
  REQUIRE(result.event == TransportEvent::GEOMETRY);
  REQUIRE_FALSE(result.temperature_field_crossing);
  p.event_cross_surface(result);
  REQUIRE(p.temperature_field_bin() == 0);
  REQUIRE(p.sqrtkT() == simulation::temperature_field.get_sqrtkT(0));
  REQUIRE(fixture.mesh->lookups == 1);
}

TEST_CASE("Rotational periodic mapping discards the incoming mesh prediction")
{
  FieldTransportFixture fixture(R"(<geometry>
    <surface id="1" type="x-plane" coeffs="0" boundary="periodic"/>
    <surface id="2" type="y-plane" coeffs="0" boundary="periodic"/>
    <surface id="3" type="z-cylinder" coeffs="0 0 2" boundary="vacuum"/>
    <cell id="1" material="void" region="1 2 -3" temperature="250"/>
  </geometry>)");
  pugi::xml_document doc;
  doc.load_string(R"(<mesh id="2"><dimension>2 2 1</dimension>
    <lower_left>0 0 -1</lower_left><upper_right>2 2 1</upper_right></mesh>)");
  fixture.mesh = make_unique<CountingMesh>(doc.document_element());
  simulation::temperature_field =
    TemperatureField(fixture.mesh.get(), {300, 400, 500, 600});
  Particle p;
  fixture.initialize(p, {0.5, 1.5, 0}, {-1, 0, 0});
  fixture.mesh->lookups = 0;
  auto result = p.event_advance();
  REQUIRE(result.event == TransportEvent::GEOMETRY);
  REQUIRE(result.next_temperature_field_bin == C_NONE);
  p.event_cross_surface(result);
  REQUIRE(p.temperature_field_bin() == 1);
  REQUIRE(p.sqrtkT() == simulation::temperature_field.get_sqrtkT(1));
  REQUIRE(fixture.mesh->lookups == 1);
}

TEST_CASE("Lattice crossings reuse the field bin, including in outer space")
{
  const double start = GENERATE(-1.25, -0.25, 0.75, 1.75);
  CAPTURE(start);
  FieldTransportFixture fixture(R"(<geometry>
    <surface id="1" type="x-plane" coeffs="-2" boundary="vacuum"/>
    <surface id="2" type="x-plane" coeffs="3" boundary="vacuum"/>
    <cell id="1" fill="10" region="1 -2"/>
    <cell id="2" universe="1" material="void" temperature="250"/>
    <cell id="3" universe="2" material="void" temperature="275"/>
    <lattice id="10"><dimension>2 1 1</dimension>
      <lower_left>0 -1 -1</lower_left><pitch>1 2 2</pitch>
      <outer>2</outer><universes>1 1</universes></lattice>
  </geometry>)");
  Particle p;
  fixture.initialize(p, {start, 0, 0});
  fixture.mesh->lookups = 0;
  auto result = p.event_advance();
  REQUIRE(result.event == TransportEvent::GEOMETRY);
  p.event_cross_surface(result);
  const int bin = start < -1 || start > 1 ? C_NONE : start < 0 ? 0 : 1;
  REQUIRE(p.temperature_field_bin() == bin);
  const double expected = bin == C_NONE
                            ? std::sqrt(275 * K_BOLTZMANN)
                            : simulation::temperature_field.get_sqrtkT(bin);
  REQUIRE(p.sqrtkT() == expected);
  REQUIRE(fixture.mesh->lookups == 0);
}
