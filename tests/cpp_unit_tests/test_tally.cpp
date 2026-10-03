#include "openmc/particle.h"
#include "openmc/tallies/filter_energy.h"
#include "openmc/tallies/filter_legendre.h"
#include "openmc/tallies/filter_lifetime_moment.h"
#include "openmc/tallies/tally.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <fmt/core.h>
#include <pugixml.hpp>

#include <stdexcept>
#include <string>

using namespace openmc;

TEST_CASE("Invalid tally nuclide names remain catchable")
{
  Tally* tally = Tally::create();
  REQUIRE_THROWS_WITH(tally->set_nuclides({"not-a-nuclide"}),
    "Could not add nuclide 'not-a-nuclide' to a tally: "
    "Nuclide 'not-a-nuclide' is not present in library.");
}

TEST_CASE("Test add/set_filter")
{
  // create a new tally object
  Tally* tally = Tally::create();

  // create a new particle filter
  Filter* particle_filter = Filter::create("particle");

  // add the particle filter to the tally
  tally->add_filter(particle_filter);

  // the filter should be added to the tally
  REQUIRE(tally->filters().size() == 1);
  REQUIRE(model::filter_map[particle_filter->id()] == tally->filters(0));

  // add the particle filter to the tally again
  tally->add_filter(particle_filter);
  // the tally should have the same number of filters
  REQUIRE(tally->filters().size() == 1);

  // create a cell filter
  Filter* cell_filter = Filter::create("cell");
  tally->add_filter(cell_filter);

  // now the size of the filters should have increased
  REQUIRE(tally->filters().size() == 2);
  REQUIRE(model::filter_map[cell_filter->id()] == tally->filters(1));

  // if we set the filters explicitly there shouldn't be extra filters hanging
  // around
  tally->set_filters({&cell_filter, 1});

  REQUIRE(tally->filters().size() == 1);
  REQUIRE(model::filter_map[cell_filter->id()] == tally->filters(0));

  // set filters again using both filters
  std::vector<Filter*> filters = {cell_filter, particle_filter};
  tally->set_filters(filters);

  REQUIRE(tally->filters().size() == 2);
  REQUIRE(model::filter_map[cell_filter->id()] == tally->filters(0));
  REQUIRE(model::filter_map[particle_filter->id()] == tally->filters(1));

  // set filters with a duplicate filter, should only add the filter to the
  // tally once
  filters = {cell_filter, cell_filter};
  tally->set_filters(filters);
  REQUIRE(tally->filters().size() == 1);
  REQUIRE(model::filter_map[cell_filter->id()] == tally->filters(0));
}

// Regression test for 64-bit tally filter-bin counts (mesh x groups > 2^31).
TEST_CASE("Tally filter-bin count does not overflow 32 bits")
{
  // Two energy filters whose bin counts multiply to 2.5e9, above INT32_MAX.
  constexpr int64_t bins_per_filter = 50000;

  // Only the bin count matters here, so the edge values are an arbitrary ramp.
  std::vector<double> edges(bins_per_filter + 1);
  for (int64_t i = 0; i < bins_per_filter + 1; ++i)
    edges[i] = static_cast<double>(i);

  Tally* tally = Tally::create();
  for (int i = 0; i < 2; ++i) {
    Filter* filter = Filter::create("energy");
    dynamic_cast<EnergyFilter*>(filter)->set_bins(edges);
    tally->add_filter(filter);
  }
  tally->set_strides();

  // set_strides() previously accumulated this product in a 32-bit int.
  REQUIRE(tally->n_filter_bins() == bins_per_filter * bins_per_filter);
  REQUIRE(tally->n_filter_bins() > 2147483647);
}

TEST_CASE("Lifetime moment filter weights are powers of the lifetime")
{
  auto* filter =
    dynamic_cast<LifetimeMomentFilter*>(Filter::create("lifetimemoment"));
  REQUIRE(filter);
  REQUIRE(filter->type_str() == "lifetimemoment");

  // A new filter holds only the zeroth moment
  REQUIRE(filter->order() == 0);
  REQUIRE(filter->n_bins() == 1);

  filter->set_order(3);
  REQUIRE(filter->order() == 3);
  REQUIRE(filter->n_bins() == 4);

  // A negative order is rejected and leaves the filter unchanged
  REQUIRE_THROWS_AS(filter->set_order(-1), std::invalid_argument);
  REQUIRE(filter->n_bins() == 4);

  // Bin n carries the weight tau^n, with tau the time since birth
  Particle p;
  const double tau = 3.0e-5;
  p.lifetime() = tau;
  FilterMatch match;
  filter->get_all_bins(p, TallyEstimator::COLLISION, match);
  REQUIRE(match.bins_ == vector<int> {0, 1, 2, 3});
  REQUIRE(
    match.weights_ == vector<double> {1.0, tau, tau * tau, tau * tau * tau});

  // At birth only the zeroth moment is nonzero
  p.lifetime() = 0.0;
  FilterMatch match_birth;
  filter->get_all_bins(p, TallyEstimator::ANALOG, match_birth);
  REQUIRE(match_birth.weights_ == vector<double> {1.0, 0.0, 0.0, 0.0});
}

TEST_CASE("Lifetime moment filter avoids the track-length estimator")
{
  auto* lifetime_filter = Filter::create("lifetimemoment");
  auto* legendre_filter =
    dynamic_cast<LegendreFilter*>(Filter::create("legendre"));
  legendre_filter->set_order(1);

  // Build a tally from XML with the given filters and estimator element
  auto tally_estimator = [](int32_t id, const std::string& filters,
                           const std::string& estimator) {
    std::string xml = fmt::format(
      "<tally id=\"{}\"><filters>{}</filters><scores>total</scores>{}</tally>",
      id, filters, estimator);
    pugi::xml_document doc;
    REQUIRE(doc.load_string(xml.c_str()));
    Tally tally(doc.child("tally"));
    return tally.estimator_;
  };
  auto lifetime_id = std::to_string(lifetime_filter->id());
  auto both_ids = fmt::format("{} {}", legendre_filter->id(), lifetime_id);

  // The default track-length estimator is replaced by a collision estimator
  REQUIRE(tally_estimator(101, lifetime_id, "") == TallyEstimator::COLLISION);

  // Explicit analog and collision estimators are honored
  REQUIRE(tally_estimator(102, lifetime_id, "<estimator>analog</estimator>") ==
          TallyEstimator::ANALOG);
  REQUIRE(tally_estimator(103, lifetime_id,
            "<estimator>collision</estimator>") == TallyEstimator::COLLISION);

  // An explicit track-length estimator is an error
  REQUIRE_THROWS_AS(
    tally_estimator(104, lifetime_id, "<estimator>tracklength</estimator>"),
    std::runtime_error);

  // An analog estimator required by another filter is kept
  REQUIRE(tally_estimator(105, both_ids, "") == TallyEstimator::ANALOG);
}

TEST_CASE("Run-time tallies with a lifetime moment filter avoid track length")
{
  auto* lifetime_filter = Filter::create("lifetimemoment");
  auto* energy_filter = Filter::create("energy");

  // A tally created at run time keeps its estimator until it is checked
  Tally* tally = Tally::create();
  tally->add_filter(lifetime_filter);
  REQUIRE(tally->estimator_ == TallyEstimator::TRACKLENGTH);
  tally->check_estimator();
  REQUIRE(tally->estimator_ == TallyEstimator::COLLISION);

  // An analog estimator is kept
  Tally* analog = Tally::create();
  analog->add_filter(lifetime_filter);
  analog->estimator_ = TallyEstimator::ANALOG;
  analog->check_estimator();
  REQUIRE(analog->estimator_ == TallyEstimator::ANALOG);

  // A mesh surface tally does not use the estimator and is not changed
  Tally* surface = Tally::create();
  surface->add_filter(lifetime_filter);
  surface->type_ = TallyType::MESH_SURFACE;
  surface->check_estimator();
  REQUIRE(surface->estimator_ == TallyEstimator::TRACKLENGTH);

  // Other tallies are not changed
  Tally* other = Tally::create();
  other->add_filter(energy_filter);
  other->check_estimator();
  REQUIRE(other->estimator_ == TallyEstimator::TRACKLENGTH);
}
