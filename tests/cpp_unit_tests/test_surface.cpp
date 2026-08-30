#include <memory>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <pugixml.hpp>

#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/surface.h"

using namespace openmc;

namespace {

template<typename T>
std::unique_ptr<T> make_surface(
  pugi::xml_document& doc, int id, const char* type, const char* coeffs)
{
  pugi::xml_node n = doc.append_child("surface");
  n.append_attribute("id") = id;
  n.append_attribute("type") = type;
  n.append_attribute("coeffs") = coeffs;
  return std::make_unique<T>(n);
}

// Register a surface under the given 1-based index so that Region can find it
template<typename T>
void add_surface(
  pugi::xml_document& doc, int id, const char* type, const char* coeffs)
{
  model::surfaces.push_back(make_surface<T>(doc, id, type, coeffs));
  model::surface_map[id] = id - 1;
}

// Builds the cell from the model attached to issue #2632
class Issue2632Fixture {
public:
  Issue2632Fixture()
  {
    // s19 (1), s48 (2), s58 (3), s62 (4), s64 (5), s65 (6), s68 (7), s101 (8)
    add_surface<SurfaceYCylinder>(doc_, 1, "y-cylinder", "0.0 0.0 17.7");
    add_surface<SurfacePlane>(doc_, 2, "plane",
      "0.7071067811865476 6.123233995736766e-17 0.7071067811865476 11.45");
    add_surface<SurfacePlane>(doc_, 3, "plane",
      "0.7071067811865476 6.123233995736766e-17 0.7071067811865476 14.35");
    add_surface<SurfacePlane>(doc_, 4, "plane",
      "6.123233995736766e-17 1.0 6.123233995736766e-17 1.5999999999999999");
    add_surface<SurfacePlane>(
      doc_, 5, "plane", "6.123233995736766e-17 1.0 6.123233995736766e-17 -1.3");
    add_surface<SurfacePlane>(doc_, 6, "plane",
      "-0.7071067811865475 6.123233995736766e-17 0.7071067811865476 1.45");
    add_surface<SurfacePlane>(doc_, 7, "plane",
      "-0.7071067811865475 6.123233995736766e-17 0.7071067811865476 -1.45");
    add_surface<SurfaceYPlane>(doc_, 8, "y-plane", "5.6");
  }

  ~Issue2632Fixture()
  {
    model::surfaces.clear();
    model::surface_map.clear();
  }

private:
  pugi::xml_document doc_;
};

} // anonymous namespace

TEST_CASE("General plane bounding box")
{
  pugi::xml_document doc;

  SECTION("Exactly axis-aligned planes bound one axis only")
  {
    // +x normal: the positive half-space starts at x = D/A
    auto px = make_surface<SurfacePlane>(doc, 1, "plane", "1.0 0.0 0.0 5.0");
    BoundingBox pos = px->bounding_box(true);
    CHECK(pos.min.x == Catch::Approx(5.0));
    CHECK(pos.min.y == -INFTY);
    CHECK(pos.min.z == -INFTY);
    CHECK(pos.max.x == INFTY);

    BoundingBox neg = px->bounding_box(false);
    CHECK(neg.max.x == Catch::Approx(5.0));
    CHECK(neg.min.x == -INFTY);
    CHECK(neg.max.y == INFTY);

    // -y normal: the sense of the bound flips with the sign of the normal
    auto ny = make_surface<SurfacePlane>(doc, 2, "plane", "0.0 -1.0 0.0 3.0");
    BoundingBox ny_pos = ny->bounding_box(true);
    CHECK(ny_pos.max.y == Catch::Approx(-3.0));
    CHECK(ny_pos.min.y == -INFTY);
    CHECK(ny_pos.min.x == -INFTY);

    BoundingBox ny_neg = ny->bounding_box(false);
    CHECK(ny_neg.min.y == Catch::Approx(-3.0));
    CHECK(ny_neg.max.y == INFTY);
  }

  SECTION("Coefficients need not be normalized")
  {
    // 4z - 10 = 0 is the same plane as z = 2.5
    auto pz = make_surface<SurfacePlane>(doc, 3, "plane", "0.0 0.0 4.0 10.0");
    CHECK(pz->bounding_box(true).min.z == Catch::Approx(2.5));
    CHECK(pz->bounding_box(false).max.z == Catch::Approx(2.5));
  }

  SECTION("Rotation roundoff is still axis aligned (issue #2632)")
  {
    // The surface s64 from the reported model: a y-plane at -1.3 written with
    // cos(pi/2) in the x and z slots.
    auto p = make_surface<SurfacePlane>(
      doc, 4, "plane", "6.123233995736766e-17 1.0 6.123233995736766e-17 -1.3");
    BoundingBox pos = p->bounding_box(true);
    CHECK(pos.min.y == Catch::Approx(-1.3));
    // The two off-axis directions must stay unbounded
    CHECK(pos.min.x == -INFTY);
    CHECK(pos.min.z == -INFTY);
    CHECK(pos.max.x == INFTY);
    CHECK(pos.max.z == INFTY);
  }

  SECTION("Oblique planes bound nothing")
  {
    auto p = make_surface<SurfacePlane>(
      doc, 5, "plane", "0.7071067811865476 0.0 0.7071067811865476 11.45");
    for (bool side : {false, true}) {
      BoundingBox bb = p->bounding_box(side);
      CHECK(bb.min.x == -INFTY);
      CHECK(bb.min.y == -INFTY);
      CHECK(bb.min.z == -INFTY);
      CHECK(bb.max.x == INFTY);
      CHECK(bb.max.y == INFTY);
      CHECK(bb.max.z == INFTY);
    }
  }

  SECTION("A plane tilted well beyond tolerance bounds nothing")
  {
    // 1e-5 is far outside PLANE_ALIGNMENT_TOL, so this must not be treated as
    // an x-plane, and in particular must not acquire a bound on y.
    auto p = make_surface<SurfacePlane>(doc, 6, "plane", "1.0 1e-5 0.0 5.0");
    BoundingBox bb = p->bounding_box(true);
    CHECK(bb.min.x == -INFTY);
    CHECK(bb.min.y == -INFTY);
  }

  SECTION("An almost-aligned plane bounds only the aligned axis")
  {
    // The off-axis coefficient is above PLANE_ALIGNMENT_TOL while the
    // normalized normal is still within it. Only x may be bounded; a spurious
    // y bound here would exclude points that lie inside the half-space.
    auto p = make_surface<SurfacePlane>(doc, 7, "plane", "1.0 1e-11 0.0 5.0");
    BoundingBox bb = p->bounding_box(true);
    CHECK(bb.min.x == Catch::Approx(5.0));
    CHECK(bb.min.y == -INFTY);
    CHECK(bb.min.z == -INFTY);
  }

  SECTION("A degenerate plane bounds nothing")
  {
    auto p = make_surface<SurfacePlane>(doc, 8, "plane", "0.0 0.0 0.0 1.0");
    BoundingBox bb = p->bounding_box(true);
    CHECK(bb.min.x == -INFTY);
    CHECK(bb.max.x == INFTY);
  }
}

TEST_CASE("Torus bounding box")
{
  pugi::xml_document doc;

  // x0 y0 z0 A B C, so the interior spans +/-B along the axis of revolution
  // and +/-(A + C) in the perpendicular directions.
  SECTION("x-torus")
  {
    auto t = make_surface<SurfaceXTorus>(
      doc, 1, "x-torus", "1.0 2.0 3.0 5.0 0.5 0.25");
    BoundingBox in = t->bounding_box(false);
    CHECK(in.min.x == Catch::Approx(0.5));
    CHECK(in.max.x == Catch::Approx(1.5));
    CHECK(in.min.y == Catch::Approx(-3.25));
    CHECK(in.max.y == Catch::Approx(7.25));
    CHECK(in.min.z == Catch::Approx(-2.25));
    CHECK(in.max.z == Catch::Approx(8.25));

    // The exterior of a torus is unbounded
    BoundingBox out = t->bounding_box(true);
    CHECK(out.min.x == -INFTY);
    CHECK(out.max.z == INFTY);
  }

  SECTION("y-torus")
  {
    auto t = make_surface<SurfaceYTorus>(
      doc, 2, "y-torus", "1.0 2.0 3.0 5.0 0.5 0.25");
    BoundingBox in = t->bounding_box(false);
    CHECK(in.min.y == Catch::Approx(1.5));
    CHECK(in.max.y == Catch::Approx(2.5));
    CHECK(in.min.x == Catch::Approx(-4.25));
    CHECK(in.max.x == Catch::Approx(6.25));
    CHECK(in.min.z == Catch::Approx(-2.25));
    CHECK(in.max.z == Catch::Approx(8.25));

    CHECK(t->bounding_box(true).max.y == INFTY);
  }

  SECTION("z-torus")
  {
    auto t = make_surface<SurfaceZTorus>(
      doc, 3, "z-torus", "1.0 2.0 3.0 5.0 0.5 0.25");
    BoundingBox in = t->bounding_box(false);
    CHECK(in.min.z == Catch::Approx(2.5));
    CHECK(in.max.z == Catch::Approx(3.5));
    CHECK(in.min.x == Catch::Approx(-4.25));
    CHECK(in.max.x == Catch::Approx(6.25));
    CHECK(in.min.y == Catch::Approx(-3.25));
    CHECK(in.max.y == Catch::Approx(7.25));

    CHECK(t->bounding_box(true).min.z == -INFTY);
  }
}

TEST_CASE("Cell bounding box matches the Python API for issue #2632")
{
  Issue2632Fixture fixture;

  // +s64 -s101 -s19 +s48 (+s58 | +s62 | -s64 | +s65 | -s68)
  Region region("5 -8 -1 2 (3 | 4 | -5 | 6 | -7)", 0);
  BoundingBox bb = region.bounding_box(0);

  // The values reported by openmc.Cell.bounding_box in the issue
  CHECK(bb.min.x == Catch::Approx(-17.7));
  CHECK(bb.min.y == Catch::Approx(-1.3));
  CHECK(bb.min.z == Catch::Approx(-17.7));
  CHECK(bb.max.x == Catch::Approx(17.7));
  CHECK(bb.max.y == Catch::Approx(5.6));
  CHECK(bb.max.z == Catch::Approx(17.7));
}
