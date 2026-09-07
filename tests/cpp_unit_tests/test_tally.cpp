#include "openmc/tallies/filter_energy.h"
#include "openmc/tallies/tally.h"
#include <catch2/catch_test_macros.hpp>

using namespace openmc;

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