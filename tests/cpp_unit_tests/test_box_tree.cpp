#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "openmc/box_tree.h"
#include "openmc/cell.h"
#include "openmc/surface.h"

#include <fmt/core.h>
#include <pugixml.hpp>

#include <algorithm>
#include <random>
#include <set>

namespace {

openmc::Position random_point(std::mt19937& rng, double lo, double hi)
{
  std::uniform_real_distribution<double> dist(lo, hi);
  return {dist(rng), dist(rng), dist(rng)};
}

openmc::Direction random_direction(std::mt19937& rng)
{
  std::normal_distribution<double> dist;
  openmc::Direction u {dist(rng), dist(rng), dist(rng)};
  return u / u.norm();
}

// Ray-box intersection by brute force
bool ray_hits(const openmc::BoundingBox& b, openmc::Position r,
  openmc::Direction u, double& t_enter)
{
  double t0 = 0.0;
  double t1 = openmc::INFTY;
  for (int i = 0; i < 3; ++i) {
    if (u[i] == 0.0) {
      if (r[i] < b.min[i] || r[i] > b.max[i])
        return false;
      continue;
    }
    double ta = (b.min[i] - r[i]) / u[i];
    double tb = (b.max[i] - r[i]) / u[i];
    t0 = std::max(t0, std::min(ta, tb));
    t1 = std::min(t1, std::max(ta, tb));
  }
  t_enter = t0;
  return t0 <= t1;
}

// Spheres with IDs 1 to n in a cube, an outer sphere with ID 1000, a z-plane
// with ID 1001 through the middle of the cube, and a distant sphere with ID
// 1002
class SpheresFixture {
public:
  explicit SpheresFixture(int n)
  {
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> radius(0.3, 1.5);
    pugi::xml_document doc;
    auto add = [&](int id, const char* type, std::string coeffs) {
      auto node = doc.append_child("surface");
      node.append_attribute("id") = id;
      node.append_attribute("type") = type;
      node.append_attribute("coeffs") = coeffs.c_str();
      openmc::model::surface_map[id] = openmc::model::surfaces.size();
      if (std::string(type) == "sphere") {
        openmc::model::surfaces.push_back(
          std::make_unique<openmc::SurfaceSphere>(node));
      } else {
        openmc::model::surfaces.push_back(
          std::make_unique<openmc::SurfaceZPlane>(node));
      }
    };
    for (int i = 1; i <= n; ++i) {
      auto c = random_point(rng, -10.0, 10.0);
      double r = radius(rng);
      centers.push_back(c);
      radii.push_back(r);
      add(i, "sphere", fmt::format("{} {} {} {}", c.x, c.y, c.z, r));
    }
    add(1000, "sphere", "0 0 0 20");
    add(1001, "z-plane", "0");
    add(1002, "sphere", "1000 0 0 1");
  }

  ~SpheresFixture()
  {
    openmc::model::surfaces.clear();
    openmc::model::surface_map.clear();
  }

  bool in_sphere(int i, openmc::Position r) const
  {
    return (r - centers[i]).norm() < radii[i];
  }

  std::vector<openmc::Position> centers;
  std::vector<double> radii;
};

} // namespace

TEST_CASE("Box tree queries match brute force")
{
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> size(0.1, 3.0);
  std::vector<openmc::BoundingBox> boxes;
  for (int i = 0; i < 500; ++i) {
    auto lo = random_point(rng, -20.0, 20.0);
    openmc::Position hi {lo.x + size(rng), lo.y + size(rng), lo.z + size(rng)};
    boxes.push_back({lo, hi});
  }
  openmc::BoxTree tree(boxes);

  SECTION("Points")
  {
    for (int k = 0; k < 2000; ++k) {
      auto r = random_point(rng, -22.0, 22.0);
      std::set<int32_t> found;
      tree.any_containing(r, [&](int32_t i) {
        found.insert(i);
        return false;
      });
      for (int i = 0; i < boxes.size(); ++i) {
        const auto& b = boxes[i];
        bool inside = r.x >= b.min.x && r.x <= b.max.x && r.y >= b.min.y &&
                      r.y <= b.max.y && r.z >= b.min.z && r.z <= b.max.z;
        // Boxes in the tree are slightly enlarged, so a box may be reported
        // for a point just outside of it, but never missed
        if (inside)
          REQUIRE(found.count(i) == 1);
      }
    }
  }

  SECTION("Rays")
  {
    for (int k = 0; k < 2000; ++k) {
      auto r = random_point(rng, -22.0, 22.0);
      auto u = random_direction(rng);
      double t_max = 15.0;
      std::set<int32_t> found;
      tree.visit_ray(r, u, t_max, [&](int32_t i, double) { found.insert(i); });
      for (int i = 0; i < boxes.size(); ++i) {
        double t;
        if (ray_hits(boxes[i], r, u, t) && t <= t_max)
          REQUIRE(found.count(i) == 1);
      }
    }
  }

  SECTION("Rays with no limit visit only the boxes they enter")
  {
    for (int k = 0; k < 2000; ++k) {
      auto r = random_point(rng, -22.0, 22.0);
      auto u = random_direction(rng);
      double t_max = openmc::INFTY;
      std::set<int32_t> found;
      tree.visit_ray(r, u, t_max, [&](int32_t i, double) { found.insert(i); });
      for (int32_t i : found) {
        // Boxes in the tree are slightly enlarged
        auto b = boxes[i];
        for (int j = 0; j < 3; ++j) {
          b.min[j] -= 1e-5;
          b.max[j] += 1e-5;
        }
        double t;
        REQUIRE(ray_hits(b, r, u, t));
      }
    }
  }

  SECTION("Axis-aligned rays")
  {
    for (int k = 0; k < 600; ++k) {
      auto r = random_point(rng, -22.0, 22.0);
      openmc::Direction u {0.0, 0.0, 0.0};
      u[k % 3] = (k / 3) % 2 == 0 ? 1.0 : -1.0;
      double t_max = openmc::INFTY;
      std::set<int32_t> found;
      tree.visit_ray(r, u, t_max, [&](int32_t i, double) { found.insert(i); });
      for (int i = 0; i < boxes.size(); ++i) {
        double t;
        if (ray_hits(boxes[i], r, u, t))
          REQUIRE(found.count(i) == 1);
      }
    }
  }
}

TEST_CASE("Regions outside of many objects")
{
  constexpr int N = 40;
  SpheresFixture fixture(N);

  // Space outside of all spheres as a simple region, and outside of the lower
  // halves of all spheres as a complex region. Both are searched with a tree.
  // Equivalent regions whose root is a union are searched without one.
  std::string simple_spec = "-1000";
  std::string complex_spec = "-1000";
  for (int i = 1; i <= N; ++i) {
    simple_spec += fmt::format(" {}", i);
    complex_spec += fmt::format(" ~(-{} -1001)", i);
  }
  // The distant sphere is not reached by the test rays
  auto untreed = [](const std::string& spec) {
    return fmt::format("({}) | -1002", spec);
  };

  struct Case {
    std::string spec;
    bool lower_halves;
  };
  for (const Case& c : {Case {simple_spec, false}, Case {complex_spec, true}}) {
    openmc::Region region(c.spec, 0);
    openmc::Region reference(untreed(c.spec), 0);
    REQUIRE_FALSE(region.is_simple());

    std::mt19937 rng(7);
    int n_inside = 0;
    for (int k = 0; k < 3000; ++k) {
      auto r = random_point(rng, -12.0, 12.0);
      auto u = random_direction(rng);

      // The region is inside sphere 1000, of radius 20
      bool expected = r.norm() < 20.0;
      for (int i = 0; i < N; ++i) {
        if (fixture.in_sphere(i, r) && (!c.lower_halves || r.z < 0.0))
          expected = false;
      }
      REQUIRE(region.contains(r, u, 0) == expected);
      if (!expected)
        continue;
      ++n_inside;

      auto [d, surf] = region.distance(r, u, 0);
      auto [d_ref, surf_ref] = reference.distance(r, u, 0);
      REQUIRE(d == Catch::Approx(d_ref).epsilon(1e-12));
      REQUIRE(surf == surf_ref);

      // Continue from the boundary: the point is outside of the region when
      // the surface is crossed in the direction of travel
      auto p = r + d * u;
      REQUIRE_FALSE(region.contains(p, u, surf));
      REQUIRE(region.contains(p, -u, -surf));
    }
    REQUIRE(n_inside > 1000);
  }
}
