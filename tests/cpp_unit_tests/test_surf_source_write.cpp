#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <pugixml.hpp>

#include "openmc/settings.h"

auto ssw_cells_from_xml(std::string xml_string)
{
  // Create a pugixml document object
  pugi::xml_document doc;

  // Load the XML from the string
  pugi::xml_parse_result result = doc.load_string(xml_string.c_str());

  pugi::xml_node root = doc.child("settings");

  // Read the surface source settings
  openmc::read_settings_xml(root);

  // Get the cell direction map
  auto ssw_cells = openmc::settings::ssw_cells;

  // Finalize to allow more tests
  openmc::free_memory_settings();
  openmc::settings::n_particles = -1;
  openmc::settings::run_mode = openmc::RunMode::UNSET;

  return ssw_cells;
}

TEST_CASE("Test duplicate cells in surface source write")
{
  std::vector<std::string> xml_strings = {
    R"(
        <settings>
          <run_mode>fixed source</run_mode>
          <particles>60</particles>
          <batches>1</batches>
          <surf_source_write>
            <surface_ids>1</surface_ids>
            <max_particles>300</max_particles>
            <cells>1 1</cells>
            <directions>to to</directions>
          </surf_source_write>
        </settings>
    )",
    R"(
        <settings>
          <run_mode>fixed source</run_mode>
          <particles>60</particles>
          <batches>1</batches>
          <surf_source_write>
            <surface_ids>1</surface_ids>
            <max_particles>300</max_particles>
            <cells>1 1</cells>
            <directions>to from</directions>
          </surf_source_write>
        </settings>
    )",
    R"(
        <settings>
          <run_mode>fixed source</run_mode>
          <particles>60</particles>
          <batches>1</batches>
          <surf_source_write>
            <surface_ids>1</surface_ids>
            <max_particles>300</max_particles>
            <cells>1 1</cells>
            <directions>to both</directions>
          </surf_source_write>
        </settings>
    )",
    R"(
        <settings>
          <run_mode>fixed source</run_mode>
          <particles>60</particles>
          <batches>1</batches>
          <surf_source_write>
            <surface_ids>1</surface_ids>
            <max_particles>300</max_particles>
            <cells>1 1</cells>
            <directions>both to</directions>
          </surf_source_write>
        </settings>
    )"};

  // to + to -> to
  auto ssw_cells = ssw_cells_from_xml(xml_strings[0]);
  REQUIRE(ssw_cells.size() == 1);
  REQUIRE(ssw_cells.at(1) == openmc::SSWCellType::To);

  // to + from -> both
  ssw_cells = ssw_cells_from_xml(xml_strings[1]);
  REQUIRE(ssw_cells.size() == 1);
  REQUIRE(ssw_cells.at(1) == openmc::SSWCellType::Both);

  // to + both -> both
  ssw_cells = ssw_cells_from_xml(xml_strings[2]);
  REQUIRE(ssw_cells.size() == 1);
  REQUIRE(ssw_cells.at(1) == openmc::SSWCellType::Both);

  // both + to -> both
  ssw_cells = ssw_cells_from_xml(xml_strings[3]);
  REQUIRE(ssw_cells.size() == 1);
  REQUIRE(ssw_cells.at(1) == openmc::SSWCellType::Both);
}
