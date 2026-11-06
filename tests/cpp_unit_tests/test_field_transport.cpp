#include <cmath>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
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

class FieldTransportFixture {
public:
  FieldTransportFixture()
  {
    settings::run_mode = RunMode::FIXED_SOURCE;
    settings::weight_window_checkpoint_surface = false;
    settings::use_shared_secondary_bank = false;
    settings::time_cutoff.fill(INFTY);
    pugi::xml_document geometry;
    geometry.load_string(R"(<geometry>
      <surface id="1" type="x-plane" coeffs="-2" boundary="vacuum"/>
      <surface id="2" type="x-plane" coeffs="2" boundary="vacuum"/>
      <cell id="1" material="void" region="1 -2" temperature="250"/>
    </geometry>)");
    read_geometry_xml(geometry.document_element());
    model::universe_level_counts.clear();
    finalize_geometry();
    finalize_cell_densities();

    pugi::xml_document doc;
    doc.load_string(R"(<mesh id="1">
      <dimension>2 1 1</dimension><lower_left>0 -1 -1</lower_left>
      <upper_right>2 1 1</upper_right></mesh>)");
    mesh = make_unique<RegularMesh>(doc.document_element());
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

  void initialize(Particle& p)
  {
    p.init_from_r_u({0.25, 0, 0}, {1, 0, 0});
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

  unique_ptr<RegularMesh> mesh;

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
