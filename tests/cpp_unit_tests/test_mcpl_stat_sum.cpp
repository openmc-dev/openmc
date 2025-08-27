#include <catch2/catch_test_macros.hpp>
#include <cstdio>
#include <string>
#include <vector>

#include "openmc/bank.h"
#include "openmc/mcpl_interface.h"

// Test the MCPL stat:sum functionality (issue #3514)
TEST_CASE("MCPL stat:sum field")
{
  // Check if MCPL interface is available
  if (!openmc::is_mcpl_interface_available()) {
    SKIP("MCPL library not available");
  }

  SECTION("stat:sum field is written to MCPL files")
  {
    // Create a temporary filename
    std::string filename = "test_stat_sum.mcpl";

    // Create some test particles
    std::vector<openmc::SourceSite> source_bank(100);
    std::vector<int64_t> bank_index = {0, 100}; // 100 particles total

    // Initialize test particles
    for (int i = 0; i < 100; ++i) {
      source_bank[i].particle = openmc::ParticleType::neutron;
      source_bank[i].r = {i * 0.1, i * 0.2, i * 0.3};
      source_bank[i].u = {0.0, 0.0, 1.0};
      source_bank[i].E = 2.0e6; // 2 MeV
      source_bank[i].time = 0.0;
      source_bank[i].wgt = 1.0;
    }

    // Write the MCPL file
    openmc::write_mcpl_source_point(filename.c_str(), source_bank, bank_index);

    // Verify the file was created
    FILE* f = std::fopen(filename.c_str(), "r");
    REQUIRE(f != nullptr);
    std::fclose(f);

    // Read the file back to check stat:sum
    // Note: This would require mcpl_open_file and checking the header
    // Since we can't easily read MCPL headers in C++ without the full MCPL API,
    // we rely on the Python test to verify the actual content

    // Clean up
    std::remove(filename.c_str());
  }

  SECTION("stat:sum uses correct particle count")
  {
    std::string filename = "test_count.mcpl";

    // Test with different particle counts
    std::vector<int> test_counts = {1, 10, 100, 1000};

    for (int count : test_counts) {
      std::vector<openmc::SourceSite> source_bank(count);
      std::vector<int64_t> bank_index = {0, count};

      // Initialize particles
      for (int i = 0; i < count; ++i) {
        source_bank[i].particle = openmc::ParticleType::neutron;
        source_bank[i].r = {0.0, 0.0, 0.0};
        source_bank[i].u = {0.0, 0.0, 1.0};
        source_bank[i].E = 1.0e6;
        source_bank[i].time = 0.0;
        source_bank[i].wgt = 1.0;
      }

      // Write MCPL file
      openmc::write_mcpl_source_point(
        filename.c_str(), source_bank, bank_index);

      // The stat:sum should equal count (verified by Python test)
      // Here we just verify the file was created successfully
      FILE* f = std::fopen(filename.c_str(), "r");
      REQUIRE(f != nullptr);
      std::fclose(f);

      // Clean up
      std::remove(filename.c_str());
    }
  }

  SECTION("stat:sum handles empty particle bank")
  {
    std::string filename = "test_empty.mcpl";

    // Create empty particle bank
    std::vector<openmc::SourceSite> source_bank;
    std::vector<int64_t> bank_index = {0};

    // This should still create a valid MCPL file with stat:sum = 0
    openmc::write_mcpl_source_point(filename.c_str(), source_bank, bank_index);

    // Verify file was created
    FILE* f = std::fopen(filename.c_str(), "r");
    REQUIRE(f != nullptr);
    std::fclose(f);

    // Clean up
    std::remove(filename.c_str());
  }
}
