#include "openmc/file_utils.h"
#include <catch2/catch_test_macros.hpp>

using namespace openmc;

TEST_CASE("Test get_file_extension")
{
#if defined(_WIN32) || defined(_WIN64)
  REQUIRE(get_file_extension("rememberthealamo.png") == "png");
  REQUIRE(get_file_extension("statepoint.20.h5") == "h5");
  REQUIRE(get_file_extension("wEiRDNaa_ame.h4") == "h4");
  REQUIRE(get_file_extension("has_directory\\asdf.20.h5") == "h5");
  REQUIRE(get_file_extension("wasssssup_lol") == "");
  REQUIRE(get_file_extension("has_directory\\secret_file") == "");
  REQUIRE(get_file_extension("lovely.dir\\extensionless_file") == "");
  REQUIRE(get_file_extension("lovely.dir\\statepoint.20.h5") == "h5");
  REQUIRE(get_file_extension("lovely.dir\\asdf123.cpp") == "cpp");

#else
  REQUIRE(get_file_extension("rememberthealamo.png") == "png");
  REQUIRE(get_file_extension("statepoint.20.h5") == "h5");
  REQUIRE(get_file_extension("wEiRDNaa_ame.h4") == "h4");
  REQUIRE(get_file_extension("has_directory/asdf.20.h5") == "h5");
  REQUIRE(get_file_extension("wasssssup_lol") == "");
  REQUIRE(get_file_extension("has_directory/secret_file") == "");
  REQUIRE(get_file_extension("lovely.dir/extensionless_file") == "");
  REQUIRE(get_file_extension("lovely.dir/statepoint.20.h5") == "h5");
  REQUIRE(get_file_extension("lovely.dir/asdf123.cpp") == "cpp");
#endif
}

TEST_CASE("Test dir_exists")
{
#if defined(_WIN32) || defined(_WIN64)
  // If this doesn't exist on a Windows system, I have no clue what is happening
  REQUIRE(dir_exists("C:\\"));

  // if this exists on your system... you deserve for this test to fail
  REQUIRE(!dir_exists("C:\\asdfa\\asdfasdf\\asdgasodgosuihasjkgh"));
#else
  // not sure how to test this when running on windows?
  REQUIRE(dir_exists("/"));

  // if this exists on your system... you deserve for this test to fail
  REQUIRE(!dir_exists("/asdfa/asdfasdf/asdgasodgosuihasjkgh/"));
#endif
}

TEST_CASE("Test file_exists")
{
#if defined(_WIN32) || defined(_WIN64)
  // Note: not clear how to portably test where a file should exist.
  REQUIRE(!file_exists("C:\\should_not_exist\\really_do_not_make_this_please"));
#else
  // Note: not clear how to portably test where a file should exist.
  REQUIRE(!file_exists("./should_not_exist/really_do_not_make_this_please"));
#endif
}

TEST_CASE("Test dir_name")
{
  REQUIRE(dir_name("") == "");
  REQUIRE(dir_name("/") == "/");
  REQUIRE(dir_name("hello") == "");
  REQUIRE(dir_name("hello/world") == "hello");
  REQUIRE(dir_name("/path/to/dir/") == "/path/to/dir");
}
