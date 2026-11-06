#include <cmath>
#include <string>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <fmt/format.h>
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
#include "openmc/tallies/tally.h"

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
    model::active_pulse_height_tallies = pulse_height_tallies_;
    model::pulse_height_cells = pulse_height_cells_;
  }

  void initialize(
    Particle& p, Position r = {0.25, 0, 0}, Direction u = {1, 0, 0})
  {
    SourceSite source;
    source.r = r;
    source.u = u;
    source.E = 1.e6;
    source.wgt = 1;
    source.particle = ParticleType::neutron();
    const int lookups = mesh->lookups;
    p.from_source(&source);
    REQUIRE(mesh->lookups == lookups);
    p.stream() = STREAM_TRACKING;
    p.seeds(0) = 123;
    p.event_calculate_xs();
    REQUIRE(p.lowest_coord().cell() != C_NONE);
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
  vector<int> pulse_height_tallies_ {model::active_pulse_height_tallies};
  vector<int32_t> pulse_height_cells_ {model::pulse_height_cells};
};

} // namespace

TEST_CASE("Advance selects the transport event and updates the field bin")
{
  FieldTransportFixture fixture;
  Particle p;
  fixture.initialize(p);
  TransportEvent expected = TransportEvent::NONE;
  int expected_bin = p.temperature_field_bin();
  double expected_x = 1;
  bool expected_alive = true;

  SECTION("Mesh crossing")
  {
    expected_bin = 1;
  }
  SECTION("Geometry without a field")
  {
    settings::temperature_field_on = false;
    expected = TransportEvent::GEOMETRY;
    expected_x = 2;
  }
  SECTION("Collision before the candidate mesh crossing")
  {
    fixture.collision_at(p, 0.25);
    expected = TransportEvent::COLLISION;
    expected_x = 0.5;
  }
  SECTION("Collision wins a tie with a mesh crossing")
  {
    fixture.collision_at(p, 0.75);
    expected = TransportEvent::COLLISION;
  }
  SECTION("Cutoff wins a tie with a mesh crossing")
  {
    settings::time_cutoff[0] = 0.75 / p.speed();
    expected_alive = false;
  }
  SECTION("Coincident geometry and field crossing")
  {
    fixture.initialize(p, {1.25, 0, 0});
    expected = TransportEvent::GEOMETRY;
    expected_bin = C_NONE;
    expected_x = 2;
  }
  SECTION("Collision wins a tie with both crossings")
  {
    fixture.initialize(p, {1.25, 0, 0});
    fixture.collision_at(p, 0.75);
    expected = TransportEvent::COLLISION;
    expected_bin = 1;
    expected_x = 2;
  }
  SECTION("Cutoff wins a tie with collision and both crossings")
  {
    fixture.initialize(p, {1.25, 0, 0});
    fixture.collision_at(p, 0.75);
    settings::time_cutoff[0] = 0.75 / p.speed();
    expected_alive = false;
    expected_bin = 1;
    expected_x = 2;
  }
  REQUIRE(p.event_advance() == expected);
  REQUIRE(p.r().x == Catch::Approx(expected_x));
  REQUIRE(p.alive() == expected_alive);
  REQUIRE(p.temperature_field_bin() == expected_bin);
}

TEST_CASE("A field-only crossing defers temperature to the XS stage")
{
  const bool event_based = GENERATE(false, true);
  const double start = GENERATE(0.25, 1.25);
  CAPTURE(event_based, start);
  FieldTransportFixture fixture(R"(<geometry>
    <surface id="1" type="x-plane" coeffs="-2" boundary="vacuum"/>
    <surface id="2" type="x-plane" coeffs="3" boundary="vacuum"/>
    <cell id="1" material="void" region="1 -2" temperature="250"/>
  </geometry>)");
  init_event_queues(1);
  auto& p = simulation::particles[0];
  fixture.initialize(p, {start, 0, 0});
  const double previous = p.sqrtkT();
  fixture.mesh->lookups = 0;
  if (event_based) {
    simulation::advance_particle_queue.thread_safe_append({p, 0});
    process_advance_particle_events();
    REQUIRE(simulation::surface_crossing_queue.size() == 0);
    REQUIRE(simulation::calculate_nonfuel_xs_queue.size() == 1);
    REQUIRE(p.n_event() == 1);
  } else {
    REQUIRE(p.event_advance() == TransportEvent::NONE);
  }
  const int bin = start < 1 ? 1 : C_NONE;
  REQUIRE(p.temperature_field_bin() == bin);
  REQUIRE(p.sqrtkT() == previous);
  if (event_based)
    process_calculate_xs_events(simulation::calculate_nonfuel_xs_queue);
  else
    p.event_calculate_xs();
  const double expected_temperature =
    bin == C_NONE ? std::sqrt(250 * K_BOLTZMANN)
                  : simulation::temperature_field.get_sqrtkT(bin);
  REQUIRE(p.sqrtkT() == expected_temperature);
  REQUIRE(p.sqrtkT_last() == previous);
  REQUIRE(fixture.mesh->lookups == 0);
}

TEST_CASE("Geometry crossings defer field temperature to the XS stage")
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
  const int bin = split < 0 || split == 2 ? C_NONE : split < 1 ? 0 : 1;
  if (event_based) {
    simulation::advance_particle_queue.thread_safe_append({p, 0});
    process_advance_particle_events();
    REQUIRE(simulation::surface_crossing_queue.size() == 1);
    REQUIRE(p.temperature_field_bin() == bin);
    REQUIRE(p.sqrtkT() == previous);
    process_surface_crossing_events();
    REQUIRE(p.n_event() == 1);
    REQUIRE(simulation::calculate_nonfuel_xs_queue.size() == 1);
  } else {
    REQUIRE(p.event_advance() == TransportEvent::GEOMETRY);
    REQUIRE(p.temperature_field_bin() == bin);
    REQUIRE(p.sqrtkT() == previous);
    p.event_cross_surface();
  }
  REQUIRE(p.lowest_coord().cell() == 1);
  REQUIRE(p.temperature_field_bin() == bin);
  REQUIRE(p.sqrtkT() == std::sqrt(275 * K_BOLTZMANN));
  REQUIRE(p.sqrtkT_last() == previous);
  if (event_based)
    process_calculate_xs_events(simulation::calculate_nonfuel_xs_queue);
  else
    p.event_calculate_xs();
  const double expected = bin == C_NONE
                            ? std::sqrt(275 * K_BOLTZMANN)
                            : simulation::temperature_field.get_sqrtkT(bin);
  REQUIRE(p.sqrtkT() == expected);
  REQUIRE(p.sqrtkT_last() == previous);
  REQUIRE(fixture.mesh->lookups == 0);
}

TEST_CASE("Boundary conditions defer field-bin correction to the XS stage")
{
  const std::string bc = GENERATE("reflective", "white", "periodic", "vacuum");
  const bool event_based = GENERATE(false, true);
  CAPTURE(bc, event_based);
  auto xml = fmt::format(R"(<geometry>
    <surface id="1" type="x-plane" coeffs="0" boundary="{}"/>
    <surface id="2" type="x-plane" coeffs="2" boundary="{}"/>
    <cell id="1" material="void" region="1 -2" temperature="250"/>
  </geometry>)",
    bc, bc);
  FieldTransportFixture fixture(xml.c_str());
  init_event_queues(1);
  auto& p = simulation::particles[0];
  fixture.initialize(p, {1.25, 0, 0});
  const double previous = p.sqrtkT();
  fixture.mesh->lookups = 0;
  if (event_based) {
    simulation::advance_particle_queue.thread_safe_append({p, 0});
    process_advance_particle_events();
    REQUIRE(simulation::surface_crossing_queue.size() == 1);
  } else {
    REQUIRE(p.event_advance() == TransportEvent::GEOMETRY);
  }
  REQUIRE(p.temperature_field_bin() == C_NONE);
  REQUIRE(p.sqrtkT() == previous);
  if (event_based)
    process_surface_crossing_events();
  else
    p.event_cross_surface();
  if (bc == "vacuum") {
    REQUIRE_FALSE(p.alive());
    REQUIRE(fixture.mesh->lookups == 0);
  } else {
    const int bin = bc == "periodic" ? 0 : 1;
    REQUIRE(fixture.mesh->lookups == 0);
    REQUIRE(p.sqrtkT() == std::sqrt(250 * K_BOLTZMANN));
    REQUIRE(p.sqrtkT_last() == previous);
    if (event_based)
      process_calculate_xs_events(simulation::calculate_nonfuel_xs_queue);
    else
      p.event_calculate_xs();
    REQUIRE(p.temperature_field_bin() == bin);
    REQUIRE(p.sqrtkT() == simulation::temperature_field.get_sqrtkT(bin));
    REQUIRE(p.sqrtkT_last() == previous);
    REQUIRE(fixture.mesh->lookups == 1);
  }
}

TEST_CASE("Periodic mapping defers relocation and caches outside results")
{
  const double lower = GENERATE(-1.0, 1.0);
  FieldTransportFixture fixture(R"(<geometry>
    <surface id="1" type="x-plane" coeffs="0" boundary="periodic"/>
    <surface id="2" type="x-plane" coeffs="2" boundary="periodic"/>
    <cell id="1" material="void" region="1 -2" temperature="250"/>
  </geometry>)");
  pugi::xml_document doc;
  auto xml = fmt::format(R"(<mesh id="2"><dimension>2 1 1</dimension>
    <lower_left>{} -1 -1</lower_left><upper_right>{} 1 1</upper_right></mesh>)",
    lower, lower + 4);
  doc.load_string(xml.c_str());
  fixture.mesh = make_unique<CountingMesh>(doc.document_element());
  simulation::temperature_field =
    TemperatureField(fixture.mesh.get(), {300, 600});
  Particle p;
  fixture.initialize(p, {1.25, 0, 0});
  const double previous = p.sqrtkT();
  fixture.mesh->lookups = 0;
  const int previous_bin = p.temperature_field_bin();
  REQUIRE(p.event_advance() == TransportEvent::GEOMETRY);
  REQUIRE(p.temperature_field_bin() == previous_bin);
  p.event_cross_surface();
  const int bin = lower < 0 ? 0 : C_NONE;
  REQUIRE(fixture.mesh->lookups == 0);
  REQUIRE(p.sqrtkT() == std::sqrt(250 * K_BOLTZMANN));
  REQUIRE(p.sqrtkT_last() == previous);
  p.event_calculate_xs();
  REQUIRE(p.temperature_field_bin() == bin);
  const double expected = bin == C_NONE
                            ? std::sqrt(250 * K_BOLTZMANN)
                            : simulation::temperature_field.get_sqrtkT(bin);
  REQUIRE(p.sqrtkT() == expected);
  REQUIRE(p.sqrtkT_last() == previous);
  REQUIRE(fixture.mesh->lookups == 1);
  p.event_calculate_xs();
  REQUIRE(p.temperature_field_bin() == bin);
  REQUIRE(p.sqrtkT() == expected);
  REQUIRE(fixture.mesh->lookups == 1);
}

TEST_CASE("Reflection at a mesh corner resolves the outgoing bin")
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
  fixture.mesh->lookups = 0;
  REQUIRE(p.event_advance() == TransportEvent::GEOMETRY);
  REQUIRE(p.temperature_field_bin() == 3);
  p.event_cross_surface();
  REQUIRE(fixture.mesh->lookups == 0);
  REQUIRE(p.sqrtkT() == std::sqrt(250 * K_BOLTZMANN));
  p.event_calculate_xs();
  REQUIRE(p.temperature_field_bin() == 2);
  REQUIRE(p.sqrtkT() == simulation::temperature_field.get_sqrtkT(2));
  REQUIRE(fixture.mesh->lookups == 1);
}

TEST_CASE("Sources initialize the field bin once, including secondary sources")
{
  const bool pulse_height = GENERATE(false, true);
  FieldTransportFixture fixture;
  if (pulse_height) {
    model::active_pulse_height_tallies = {0};
    model::pulse_height_cells = {0};
  }
  Particle p;
  fixture.initialize(p);
  REQUIRE(p.temperature_field_bin() == 0);
  REQUIRE(p.sqrtkT() == simulation::temperature_field.get_sqrtkT(0));
  REQUIRE(fixture.mesh->lookups == 1);
  fixture.mesh->lookups = 0;
  SourceSite secondary;
  secondary.r = {1.25, 0, 0};
  secondary.u = {1, 0, 0};
  secondary.E = 1.e6;
  secondary.wgt = 1;
  secondary.particle = ParticleType::photon();
  p.event_revive_from_secondary(secondary);
  REQUIRE(fixture.mesh->lookups == 0);
  if (pulse_height) {
    REQUIRE(p.lowest_coord().cell() == 0);
    REQUIRE(p.sqrtkT() == std::sqrt(250 * K_BOLTZMANN));
    REQUIRE(p.pht_storage()[0] == -secondary.E);
  }
  p.event_calculate_xs();
  REQUIRE(p.temperature_field_bin() == 1);
  REQUIRE(p.sqrtkT() == simulation::temperature_field.get_sqrtkT(1));
  REQUIRE(fixture.mesh->lookups == 1);
}

TEST_CASE("Geometry queries preserve the bin and assign cell temperatures")
{
  FieldTransportFixture fixture;
  GeometryState p;
  p.init_from_r_u({1.25, 0, 0}, {1, 0, 0});
  p.temperature_field_bin() = 0;
  REQUIRE(exhaustive_find_cell(p));
  REQUIRE(p.temperature_field_bin() == 0);
  REQUIRE(p.sqrtkT() == std::sqrt(250 * K_BOLTZMANN));
  REQUIRE(fixture.mesh->lookups == 0);
}

TEST_CASE("Lattice crossings defer field temperature to the XS stage")
{
  const double start = GENERATE(-1.25, -0.25, 0.75, 1.75);
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
  REQUIRE(p.event_advance() == TransportEvent::GEOMETRY);
  p.event_cross_surface();
  const int bin = start < -1 || start > 1 ? C_NONE : start < 0 ? 0 : 1;
  REQUIRE(p.temperature_field_bin() == bin);
  p.event_calculate_xs();
  const double expected = bin == C_NONE
                            ? std::sqrt(275 * K_BOLTZMANN)
                            : simulation::temperature_field.get_sqrtkT(bin);
  REQUIRE(p.sqrtkT() == expected);
  REQUIRE(fixture.mesh->lookups == 0);
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
