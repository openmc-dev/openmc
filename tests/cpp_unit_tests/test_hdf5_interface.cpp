#include <catch2/catch_test_macros.hpp>

#include <cstdio>
#include <string>

#include "openmc/hdf5_interface.h"
#include "openmc/vector.h"

using namespace openmc;

namespace {

//! Scoped HDF5 file that is created on construction and deleted on destruction
class TempFile {
public:
  explicit TempFile(const std::string& filename) : filename_(filename)
  {
    file_id_ = file_open(filename_, 'w');
  }

  ~TempFile()
  {
    if (file_id_ >= 0)
      file_close(file_id_);
    std::remove(filename_.c_str());
  }

  //! Close the file and reopen it for reading, so that tests exercise what
  //! actually landed on disk rather than a cached in-memory value
  hid_t reopen()
  {
    file_close(file_id_);
    file_id_ = file_open(filename_, 'r');
    return file_id_;
  }

  hid_t id() const { return file_id_; }

private:
  std::string filename_;
  hid_t file_id_;
};

} // namespace

TEST_CASE("String datasets round-trip")
{
  TempFile file("test_hdf5_string_dataset.h5");

  const std::string empty {""};
  const std::string nonempty {"pincell.exo"};
  // A string whose length is exactly the datatype size, i.e. with no room for a
  // trailing null character. Reading must not truncate it.
  const std::string exact {"abc"};

  write_dataset(file.id(), "empty", empty);
  write_dataset(file.id(), "nonempty", nonempty);
  write_dataset(file.id(), "exact", exact);
  write_dataset(file.id(), "empty_literal", "");

  hid_t file_id = file.reopen();

  SECTION("An empty string is still written as a dataset")
  {
    // Regression test for openmc-dev/openmc#2285: writing nothing at all left
    // the dataset absent from the file and broke readers downstream
    REQUIRE(object_exists(file_id, "empty"));
    REQUIRE(object_exists(file_id, "empty_literal"));

    // Stored as a single null byte, since HDF5 has no zero-size string type
    REQUIRE(dataset_typesize(file_id, "empty") == 1);
  }

  SECTION("An empty string reads back as empty")
  {
    std::string value {"not empty"};
    read_dataset(file_id, "empty", value);
    REQUIRE(value.empty());
    REQUIRE(value == "");

    value = "not empty";
    read_dataset(file_id, "empty_literal", value);
    REQUIRE(value.empty());
  }

  SECTION("A non-empty string is unchanged")
  {
    std::string value;
    read_dataset(file_id, "nonempty", value);
    REQUIRE(value == nonempty);
    REQUIRE(dataset_typesize(file_id, "nonempty") == nonempty.size());
  }

  SECTION("A string with no room for a null terminator is not truncated")
  {
    std::string value;
    read_dataset(file_id, "exact", value);
    REQUIRE(value == exact);
    REQUIRE(value.size() == 3);
  }
}

TEST_CASE("String datasets round-trip within a group")
{
  TempFile file("test_hdf5_string_group.h5");

  hid_t group = create_group(file.id(), "geometry");
  write_dataset(group, "name", std::string {""});
  write_dataset(group, "region", std::string {"1 -2 3"});
  close_group(group);

  hid_t file_id = file.reopen();
  group = open_group(file_id, "geometry");

  REQUIRE(object_exists(group, "name"));

  std::string name {"stale"};
  read_dataset(group, "name", name);
  REQUIRE(name.empty());

  std::string region;
  read_dataset(group, "region", region);
  REQUIRE(region == "1 -2 3");

  close_group(group);
}

TEST_CASE("String attributes round-trip")
{
  TempFile file("test_hdf5_string_attribute.h5");

  write_attribute(file.id(), "empty", std::string {""});
  write_attribute(file.id(), "nonempty", std::string {"/path/to/inputs/"});
  write_attribute(file.id(), "empty_literal", "");
  write_attribute(file.id(), "exact", std::string {"abc"});

  hid_t file_id = file.reopen();

  SECTION("An empty attribute is still written")
  {
    // settings::path_input is empty unless -i is passed, so this is the common
    // case for the "path" attribute of a statepoint file
    REQUIRE(attribute_exists(file_id, "empty"));
    REQUIRE(attribute_exists(file_id, "empty_literal"));
    REQUIRE(attribute_typesize(file_id, "empty") == 1);
  }

  SECTION("An empty attribute reads back as empty")
  {
    std::string value {"not empty"};
    read_attribute(file_id, "empty", value);
    REQUIRE(value.empty());

    value = "not empty";
    read_attribute(file_id, "empty_literal", value);
    REQUIRE(value.empty());
  }

  SECTION("A non-empty attribute is unchanged")
  {
    std::string value;
    read_attribute(file_id, "nonempty", value);
    REQUIRE(value == "/path/to/inputs/");

    read_attribute(file_id, "exact", value);
    REQUIRE(value == "abc");
  }
}

TEST_CASE("Vectors of strings round-trip")
{
  TempFile file("test_hdf5_string_vector.h5");

  const vector<std::string> mixed {"U235", "", "H1"};
  const vector<std::string> all_empty {"", "", ""};
  const vector<std::string> none {};

  write_dataset(file.id(), "mixed", mixed);
  write_dataset(file.id(), "all_empty", all_empty);
  write_dataset(file.id(), "none", none);

  hid_t file_id = file.reopen();

  SECTION("Individual empty entries do not shrink the datatype")
  {
    REQUIRE(object_exists(file_id, "mixed"));
    // Sized by the longest entry plus a null terminator
    REQUIRE(dataset_typesize(file_id, "mixed") == 5);

    hid_t dset = open_dataset(file_id, "mixed");
    REQUIRE(object_shape(dset)[0] == 3);
    close_dataset(dset);
  }

  SECTION("A vector of empty strings is written")
  {
    REQUIRE(object_exists(file_id, "all_empty"));
    REQUIRE(dataset_typesize(file_id, "all_empty") == 1);

    hid_t dset = open_dataset(file_id, "all_empty");
    REQUIRE(object_shape(dset)[0] == 3);
    close_dataset(dset);
  }

  SECTION("An empty vector is written as a zero-length dataset")
  {
    REQUIRE(object_exists(file_id, "none"));

    hid_t dset = open_dataset(file_id, "none");
    REQUIRE(object_shape(dset)[0] == 0);
    close_dataset(dset);
  }
}
