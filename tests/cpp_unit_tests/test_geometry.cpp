#include <catch2/catch_test_macros.hpp>

#include "openmc/cell.h"
#include "openmc/geometry.h"
#include "openmc/particle_data.h"
#include "openmc/surface.h"
#include "openmc/universe.h"

#include <memory>
#include <string>

#include <pugixml.hpp>

namespace {

class GeometryFixture {
public:
  GeometryFixture()
    : root_universe_ {openmc::model::root_universe},
      n_coord_levels_ {openmc::model::n_coord_levels}
  {
    openmc::model::cells.clear();
    openmc::model::cell_map.clear();
    openmc::model::universes.clear();
    openmc::model::universe_map.clear();
    openmc::model::surfaces.clear();
    openmc::model::surface_map.clear();

    openmc::model::n_coord_levels = 2;
    openmc::model::root_universe = 0;

    pugi::xml_document surface_doc;
    auto surface_node = surface_doc.append_child("surface");
    surface_node.append_attribute("id") = 1;
    surface_node.append_attribute("type") = "sphere";
    surface_node.append_attribute("coeffs") = "0 0 0 1";
    openmc::model::surfaces.push_back(
      std::make_unique<openmc::SurfaceSphere>(surface_node));
    openmc::model::surface_map[1] = 0;

    openmc::model::cells.push_back(make_cell(1, 0, 1, ""));
    openmc::model::cells.push_back(make_cell(2, 1, -1, "-1"));
    openmc::model::cells.push_back(make_cell(3, 1, -1, "+1"));
    for (int i = 0; i < openmc::model::cells.size(); ++i)
      openmc::model::cell_map[openmc::model::cells[i]->id_] = i;

    auto root = std::make_unique<openmc::Universe>();
    root->id_ = 0;
    root->cells_ = {0};
    root->n_instances_ = 1;
    openmc::model::universes.push_back(std::move(root));
    openmc::model::universe_map[0] = 0;

    auto nested = std::make_unique<openmc::Universe>();
    nested->id_ = 1;
    nested->cells_ = {1, 2};
    nested->n_instances_ = 1;
    openmc::model::universes.push_back(std::move(nested));
    openmc::model::universe_map[1] = 1;
  }

  ~GeometryFixture()
  {
    openmc::model::cells.clear();
    openmc::model::cell_map.clear();
    openmc::model::universes.clear();
    openmc::model::universe_map.clear();
    openmc::model::surfaces.clear();
    openmc::model::surface_map.clear();
    openmc::model::root_universe = root_universe_;
    openmc::model::n_coord_levels = n_coord_levels_;
  }

private:
  static std::unique_ptr<openmc::CSGCell> make_cell(
    int id, int universe, int fill, const char* region)
  {
    pugi::xml_document doc;
    auto node = doc.append_child("cell");
    node.append_attribute("id") = id;
    node.append_attribute("universe") = universe;
    if (fill >= 0) {
      const auto fill_value {std::to_string(fill)};
      node.append_child("fill").text() = fill_value.c_str();
    } else {
      node.append_child("material").text() = "void";
    }
    if (region[0] != '\0')
      node.append_child("region").text() = region;

    auto cell = std::make_unique<openmc::CSGCell>(node);
    if (fill < 0) {
      cell->type_ = openmc::Fill::MATERIAL;
      cell->sqrtkT_.push_back(0.0);
      cell->density_mult_.push_back(1.0);
    } else {
      cell->type_ = openmc::Fill::UNIVERSE;
    }
    return cell;
  }

  int root_universe_;
  int n_coord_levels_;
};

} // namespace

TEST_CASE("Reconcile a particle after a collision near a surface")
{
  GeometryFixture fixture;

  openmc::GeometryState p;
  p.n_coord() = 2;
  p.coord(0).universe() = 0;
  p.coord(0).cell() = 0;
  p.coord(0).r() = {0.0, 0.0, 0.0};
  p.coord(0).u() = {-1.0, 0.0, 0.0};
  p.coord(1).universe() = 1;
  p.coord(1).cell() = 2;
  p.coord(1).r() = {1.0 - 1.0e-13, 0.0, 0.0};
  p.coord(1).u() = {-1.0, 0.0, 0.0};

  openmc::reconcile_cell_after_collision(p);
  REQUIRE(p.n_coord() == 2);
  REQUIRE(p.coord(0).cell() == 0);
  REQUIRE(p.coord(1).cell() == 1);
}
