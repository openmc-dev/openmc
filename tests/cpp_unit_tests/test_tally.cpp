#include "openmc/capi.h"
#include "openmc/particle.h"
#include "openmc/tallies/filter_energy.h"
#include "openmc/tallies/filter_legendre.h"
#include "openmc/tallies/filter_lifetime_moment.h"
#include "openmc/tallies/filter_sph_harm.h"
#include "openmc/tallies/tally.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <fmt/core.h>
#include <pugixml.hpp>

#include <iostream>
#include <sstream>
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

namespace {

// Redirects std::cerr, to which warning() writes, while the object exists
class CerrCapture {
public:
  CerrCapture() : saved_ {std::cerr.rdbuf(buffer_.rdbuf())} {}
  ~CerrCapture() { std::cerr.rdbuf(saved_); }

  // The captured text, with each run of white space (including the line
  // breaks and indentation added by warning()) replaced by one space
  std::string text() const
  {
    std::istringstream stream {buffer_.str()};
    std::string word;
    std::string result;
    while (stream >> word)
      result += (result.empty() ? "" : " ") + word;
    return result;
  }

private:
  std::ostringstream buffer_;
  std::streambuf* saved_;
};

// Check the estimator of a tally as simulation initialization does, and
// return the warnings written
std::string check_estimator(Tally* tally)
{
  CerrCapture capture;
  tally->check_estimator();
  return capture.text();
}

// Number of occurrences of a string in a text
int count(const std::string& text, const std::string& sub)
{
  int n = 0;
  for (auto pos = text.find(sub); pos != std::string::npos;
       pos = text.find(sub, pos + sub.size()))
    ++n;
  return n;
}

// One warning that names the tally and the estimator that it now uses
void check_warning(const std::string& text, const Tally& tally)
{
  bool analog = tally.estimator_ == TallyEstimator::ANALOG;
  CHECK(count(text, "WARNING:") == 1);
  CHECK(count(text, fmt::format("Tally {} has a lifetime moment filter",
                      tally.id())) == 1);
  CHECK(count(text, "An analog estimator") == (analog ? 1 : 0));
  CHECK(count(text, "A collision estimator") == (analog ? 0 : 1));
}

std::string filter_types(const vector<Filter*>& filters)
{
  std::string types;
  for (const auto* f : filters)
    types += " " + f->type_str();
  return types;
}

} // namespace

TEST_CASE("Run-time lifetime moment tallies get the estimator of tallies.xml")
{
  auto* lifetime = Filter::create("lifetimemoment");
  auto* energyout = Filter::create("energyout");
  auto* energy = Filter::create("energy");
  auto* legendre = dynamic_cast<LegendreFilter*>(Filter::create("legendre"));
  legendre->set_order(1);
  auto* production = Filter::create("particleproduction");
  auto* spatial = Filter::create("spatiallegendre");
  auto* harmonics_scatter = dynamic_cast<SphericalHarmonicsFilter*>(
    Filter::create("sphericalharmonics"));
  harmonics_scatter->set_order(1);
  harmonics_scatter->set_cosine("scatter");
  auto* harmonics_particle = dynamic_cast<SphericalHarmonicsFilter*>(
    Filter::create("sphericalharmonics"));
  harmonics_particle->set_order(1);
  harmonics_particle->set_cosine("particle");

  using TE = TallyEstimator;
  struct Case {
    vector<Filter*> filters;
    std::string score;
    TallyEstimator estimator;   // of the same tally read from tallies.xml
    TallyEstimator from_scores; // chosen by set_scores() after the filters
  };
  vector<Case> cases {
    // Another filter requires an analog estimator
    {{energyout, lifetime}, "nu-fission", TE::ANALOG, TE::ANALOG},
    {{lifetime, energyout}, "nu-fission", TE::ANALOG, TE::ANALOG},
    {{legendre, lifetime}, "scatter", TE::ANALOG, TE::ANALOG},
    {{legendre, lifetime}, "flux", TE::ANALOG, TE::TRACKLENGTH},
    {{lifetime, legendre}, "flux", TE::ANALOG, TE::TRACKLENGTH},
    {{lifetime, legendre, energy}, "flux", TE::ANALOG, TE::TRACKLENGTH},
    {{production, lifetime}, "events", TE::ANALOG, TE::TRACKLENGTH},
    {{lifetime, production}, "events", TE::ANALOG, TE::TRACKLENGTH},
    {{harmonics_scatter, lifetime}, "scatter", TE::ANALOG, TE::TRACKLENGTH},
    {{lifetime, harmonics_scatter}, "scatter", TE::ANALOG, TE::TRACKLENGTH},
    // No other filter requires an analog estimator
    {{harmonics_particle, lifetime}, "scatter", TE::COLLISION, TE::TRACKLENGTH},
    {{lifetime, harmonics_particle}, "scatter", TE::COLLISION, TE::TRACKLENGTH},
    {{spatial, lifetime}, "flux", TE::COLLISION, TE::TRACKLENGTH},
    {{lifetime, spatial}, "flux", TE::COLLISION, TE::TRACKLENGTH},
    {{energy, lifetime}, "flux", TE::COLLISION, TE::TRACKLENGTH},
    {{lifetime}, "flux", TE::COLLISION, TE::TRACKLENGTH},
  };

  int32_t id = 1000;
  for (auto& c : cases) {
    INFO("filters:" << filter_types(c.filters) << ", score: " << c.score);

    // The same tally read from tallies.xml
    {
      std::string ids;
      for (const auto* f : c.filters)
        ids += fmt::format(" {}", f->id());
      std::string xml = fmt::format(
        "<tally id=\"{}\"><filters>{}</filters><scores>{}</scores></tally>",
        ++id, ids, c.score);
      pugi::xml_document doc;
      REQUIRE(doc.load_string(xml.c_str()));
      Tally tally(doc.child("tally"));
      CHECK(tally.estimator_ == c.estimator);
    }

    // Filters set before the scores, so that set_scores() may already choose
    // an analog estimator
    Tally* filters_first = Tally::create(++id);
    filters_first->set_filters(c.filters);
    filters_first->set_scores({c.score});
    CHECK(filters_first->estimator_ == c.from_scores);
    auto text = check_estimator(filters_first);
    CHECK(filters_first->estimator_ == c.estimator);
    if (c.from_scores == TE::TRACKLENGTH) {
      check_warning(text, *filters_first);
    } else {
      CHECK(text.empty());
    }

    // Scores set before the filters, so that only check_estimator() sees the
    // filters
    Tally* scores_first = Tally::create(++id);
    scores_first->set_scores({c.score});
    scores_first->set_filters(c.filters);
    CHECK(scores_first->estimator_ == TE::TRACKLENGTH);
    text = check_estimator(scores_first);
    CHECK(scores_first->estimator_ == c.estimator);
    check_warning(text, *scores_first);

    // Checking again changes nothing
    CHECK(check_estimator(scores_first).empty());
    CHECK(scores_first->estimator_ == c.estimator);
  }
}

TEST_CASE("Run-time lifetime moment tallies keep estimators other than track "
          "length")
{
  auto* lifetime = Filter::create("lifetimemoment");
  auto* energyout = Filter::create("energyout");
  auto* legendre = dynamic_cast<LegendreFilter*>(Filter::create("legendre"));
  legendre->set_order(1);
  auto* harmonics = dynamic_cast<SphericalHarmonicsFilter*>(
    Filter::create("sphericalharmonics"));
  harmonics->set_order(1);
  harmonics->set_cosine("particle");

  // A tally with its estimator set through the C API, as by openmc.lib
  auto create = [](int32_t id, vector<Filter*> filters,
                  const std::string& score, const char* estimator) {
    Tally* tally = Tally::create(id);
    tally->set_scores({score});
    tally->set_filters(filters);
    REQUIRE(
      openmc_tally_set_estimator(model::tally_map.at(id), estimator) == 0);
    return tally;
  };

  // Only a track-length estimator is replaced: a collision estimator set at
  // run time is kept even beside a filter that requires an analog estimator
  Tally* collision =
    create(2001, {energyout, lifetime}, "nu-fission", "collision");
  CHECK(check_estimator(collision).empty());
  CHECK(collision->estimator_ == TallyEstimator::COLLISION);

  // An analog estimator is kept, even where a collision estimator would do
  Tally* analog = create(2002, {harmonics, lifetime}, "scatter", "analog");
  CHECK(check_estimator(analog).empty());
  CHECK(analog->estimator_ == TallyEstimator::ANALOG);

  // A track-length estimator set explicitly is replaced like the default one
  Tally* tracklength =
    create(2003, {legendre, lifetime}, "flux", "tracklength");
  auto text = check_estimator(tracklength);
  CHECK(tracklength->estimator_ == TallyEstimator::ANALOG);
  check_warning(text, *tracklength);
}

TEST_CASE("Run-time tallies without a lifetime moment filter keep their "
          "estimator")
{
  auto* energyout = Filter::create("energyout");
  auto* energy = Filter::create("energy");
  auto* legendre = dynamic_cast<LegendreFilter*>(Filter::create("legendre"));
  legendre->set_order(1);
  auto* production = Filter::create("particleproduction");
  auto* spatial = Filter::create("spatiallegendre");
  auto* harmonics = dynamic_cast<SphericalHarmonicsFilter*>(
    Filter::create("sphericalharmonics"));
  harmonics->set_order(1);
  harmonics->set_cosine("scatter");

  // check_estimator() only changes tallies with a lifetime moment filter;
  // other run-time tallies keep the estimator that set_scores() chose
  using TE = TallyEstimator;
  struct Case {
    vector<Filter*> filters;
    std::string score;
    TallyEstimator from_scores; // chosen by set_scores() after the filters
  };
  vector<Case> cases {
    {{energyout}, "nu-fission", TE::ANALOG},
    {{legendre}, "flux", TE::TRACKLENGTH},
    {{production}, "events", TE::TRACKLENGTH},
    {{harmonics}, "scatter", TE::TRACKLENGTH},
    {{spatial}, "flux", TE::TRACKLENGTH},
    {{energy}, "flux", TE::TRACKLENGTH},
  };

  int32_t id = 3000;
  for (auto& c : cases) {
    INFO("filters:" << filter_types(c.filters) << ", score: " << c.score);

    Tally* filters_first = Tally::create(++id);
    filters_first->set_filters(c.filters);
    filters_first->set_scores({c.score});
    CHECK(filters_first->estimator_ == c.from_scores);
    CHECK(check_estimator(filters_first).empty());
    CHECK(filters_first->estimator_ == c.from_scores);

    Tally* scores_first = Tally::create(++id);
    scores_first->set_scores({c.score});
    scores_first->set_filters(c.filters);
    CHECK(scores_first->estimator_ == TE::TRACKLENGTH);
    CHECK(check_estimator(scores_first).empty());
    CHECK(scores_first->estimator_ == TE::TRACKLENGTH);
  }
}
