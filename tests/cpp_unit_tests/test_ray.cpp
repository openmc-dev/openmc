#include <cmath>
#include <utility>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <pugixml.hpp>

#include "openmc/geometry.h"
#include "openmc/geometry_aux.h"
#include "openmc/ray.h"
#include "openmc/settings.h"
#include "openmc/surface.h"

namespace {

using Catch::Matchers::WithinAbs;

constexpr double DISTANCE_TOLERANCE = 1.0e-12;

class SphericalShellFixture {
public:
  SphericalShellFixture()
    : run_mode_(openmc::settings::run_mode),
      root_universe_(openmc::model::root_universe),
      n_coord_levels_(openmc::model::n_coord_levels)
  {
    openmc::settings::run_mode = openmc::RunMode::PLOTTING;

    constexpr auto geometry = R"(
      <geometry>
        <surface id="1" type="sphere" coeffs="0 0 0 1"/>
        <surface id="2" type="sphere" coeffs="0 0 0 2"
                 boundary="vacuum"/>
        <cell id="1" material="void" region="-1"/>
        <cell id="2" material="void" region="1 -2"/>
      </geometry>
    )";

    auto result = document_.load_string(geometry);
    REQUIRE(result);
    openmc::read_geometry_xml(document_.document_element());
    openmc::finalize_geometry();
    openmc::finalize_cell_densities();
  }

  ~SphericalShellFixture()
  {
    openmc::free_memory_geometry();
    openmc::free_memory_surfaces();
    openmc::model::universe_level_counts.clear();
    openmc::model::root_universe = root_universe_;
    openmc::model::n_coord_levels = n_coord_levels_;
    openmc::settings::run_mode = run_mode_;
  }

private:
  pugi::xml_document document_;
  openmc::RunMode run_mode_;
  int root_universe_;
  int n_coord_levels_;
};

class RecordingRay : public openmc::Ray {
public:
  using Ray::Ray;

  void on_intersection() override
  {
    intersections.emplace_back(
      traversal_distance_, boundary().surface_index() + 1);
  }

  openmc::vector<std::pair<double, int>> intersections;
};

} // namespace

TEST_CASE_METHOD(SphericalShellFixture, "Trace a single chord through a sphere")
{
  constexpr double impact_parameter = 1.5;
  const double half_chord =
    std::sqrt(4.0 - impact_parameter * impact_parameter);

  RecordingRay ray({-3.0, impact_parameter, 0.0}, {1.0, 0.0, 0.0});
  ray.trace();

  REQUIRE(ray.intersections.size() == 2);
  CHECK(ray.intersections[0].second == 2);
  CHECK_THAT(ray.intersections[0].first, WithinAbs(0.0, DISTANCE_TOLERANCE));
  CHECK(ray.intersections[1].second == 2);
  CHECK_THAT(ray.intersections[1].first,
    WithinAbs(2.0 * half_chord, DISTANCE_TOLERANCE));
}

TEST_CASE_METHOD(
  SphericalShellFixture, "Trace both segments of a spherical shell")
{
  constexpr double impact_parameter = 0.5;
  const double outer = std::sqrt(4.0 - impact_parameter * impact_parameter);
  const double inner = std::sqrt(1.0 - impact_parameter * impact_parameter);

  RecordingRay ray({-3.0, impact_parameter, 0.0}, {1.0, 0.0, 0.0});
  ray.trace();

  REQUIRE(ray.intersections.size() == 4);
  CHECK(ray.intersections[0].second == 2);
  CHECK_THAT(ray.intersections[0].first, WithinAbs(0.0, DISTANCE_TOLERANCE));
  CHECK(ray.intersections[1].second == 1);
  CHECK_THAT(
    ray.intersections[1].first, WithinAbs(outer - inner, DISTANCE_TOLERANCE));
  CHECK(ray.intersections[2].second == 1);
  CHECK_THAT(
    ray.intersections[2].first, WithinAbs(outer + inner, DISTANCE_TOLERANCE));
  CHECK(ray.intersections[3].second == 2);
  CHECK_THAT(
    ray.intersections[3].first, WithinAbs(2.0 * outer, DISTANCE_TOLERANCE));
}
