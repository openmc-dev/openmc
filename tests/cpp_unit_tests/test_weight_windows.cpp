#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <pugixml.hpp>

#include "openmc/geometry.h"
#include "openmc/mesh.h"
#include "openmc/particle.h"
#include "openmc/settings.h"
#include "openmc/weight_windows.h"

using namespace openmc;

TEST_CASE("Weight window reference is set by the first valid window")
{
  const bool weight_windows_on = settings::weight_windows_on;
  const int max_history_splits = settings::max_history_splits;
  const int n_coord_levels = model::n_coord_levels;

  struct Cleanup {
    bool weight_windows_on;
    int max_history_splits;
    int n_coord_levels;

    ~Cleanup()
    {
      free_memory_weight_windows();
      model::mesh_map.clear();
      model::meshes.clear();
      settings::weight_windows_on = weight_windows_on;
      settings::max_history_splits = max_history_splits;
      model::n_coord_levels = n_coord_levels;
    }
  } cleanup {weight_windows_on, max_history_splits, n_coord_levels};

  settings::weight_windows_on = true;
  settings::max_history_splits = 10;
  model::n_coord_levels = 1;

  pugi::xml_document mesh_doc;
  auto mesh_node = mesh_doc.append_child("mesh");
  mesh_node.append_attribute("id") = 1;
  mesh_node.append_child("dimension").text() = "5 1 1";
  mesh_node.append_child("lower_left").text() = "0 -1 -1";
  mesh_node.append_child("upper_right").text() = "5 1 1";
  model::meshes.push_back(make_unique<RegularMesh>(mesh_node));
  model::mesh_map[1] = 0;

  auto* ww = WeightWindows::create();
  vector<double> energy_bounds {0.0, 20.0e6};
  ww->set_energy_bounds(energy_bounds);
  ww->set_mesh(model::meshes.front());
  vector<double> lower_bounds {-1.0, 0.0, 1.0e-6, 0.0, 1.0e-6};
  vector<double> upper_bounds {-1.0, 0.0, 5.0e-6, 0.0, 2.0e-6};
  ww->set_bounds(lower_bounds, upper_bounds);

  auto check_first_valid_window = [ww](double first_x, ParticleType type) {
    ww->set_particle_type(type);
    Particle p;
    p.type() = type;
    p.E() = 1.0e6;
    p.wgt() = 2.0;
    p.wgt_born() = 2.0;
    p.r() = {first_x, 0.0, 0.0};

    // A matching invalid bound, or no mesh match, must leave the reference
    // unset until this particle reaches a valid window.
    apply_weight_windows(p);
    CHECK(p.wgt_ww_born() == -1.0);
    CHECK(p.wgt() == Catch::Approx(2.0));
    CHECK(p.n_split() == 0);

    // The first valid window is tiny, but normalization by its midpoint keeps
    // this non-unit birth weight inside the window instead of splitting it.
    p.r() = {2.5, 0.0, 0.0};
    apply_weight_windows(p);
    CHECK(p.wgt_ww_born() == Catch::Approx(3.0e-6));
    CHECK(p.wgt() == Catch::Approx(2.0));
    CHECK(p.n_split() == 0);
    CHECK(p.local_secondary_bank().empty());

    // A later invalid window must preserve the first valid reference.
    p.r() = {3.5, 0.0, 0.0};
    apply_weight_windows(p);
    CHECK(p.wgt_ww_born() == Catch::Approx(3.0e-6));
    CHECK(p.n_split() == 0);

    // A genuinely tighter later window still requires splitting against that
    // same reference value, conserving the parent's weight across both sites.
    p.r() = {4.5, 0.0, 0.0};
    apply_weight_windows(p);
    CHECK(p.n_split() == 2);
    CHECK(p.local_secondary_bank().size() == 1);
    if (!p.local_secondary_bank().empty()) {
      CHECK(
        p.local_secondary_bank().front().wgt_ww_born == Catch::Approx(3.0e-6));
      CHECK(
        p.wgt() + p.local_secondary_bank().front().wgt == Catch::Approx(2.0));
    }
  };

  SECTION("neutron weight windows")
  {
    SECTION("matching negative lower bound")
    {
      check_first_valid_window(0.5, ParticleType::neutron());
    }

    SECTION("matching zero lower bound")
    {
      check_first_valid_window(1.5, ParticleType::neutron());
    }

    SECTION("no initial mesh match")
    {
      check_first_valid_window(-0.5, ParticleType::neutron());
    }
  }

  SECTION("photon weight windows")
  {
    SECTION("matching negative lower bound")
    {
      check_first_valid_window(0.5, ParticleType::photon());
    }

    SECTION("matching zero lower bound")
    {
      check_first_valid_window(1.5, ParticleType::photon());
    }

    SECTION("no initial mesh match")
    {
      check_first_valid_window(-0.5, ParticleType::photon());
    }
  }
}
