#ifndef OPENMC_CROSS_SECTIONS_H
#define OPENMC_CROSS_SECTIONS_H

#include "pugixml.hpp"

#include <string>
#include <map>
#include <vector>

namespace openmc {

//==============================================================================
// Library class
//==============================================================================

class Library {
public:
  // Types, enums
  enum class Type {
    neutron = 1, photon = 3, thermal = 2, multigroup = 4, wmp = 5
  };

  // Constructors
  Library() { };
  Library(pugi::xml_node node, const std::string& directory);

  // Comparison operator (for using in map)
  bool operator<(const Library& other) {
    return path_ < other.path_;
  }

  // Data members
  Type type_;
  std::vector<std::string> materials_;
  std::string path_;
};

//==============================================================================
// Global variable declarations
//==============================================================================

extern std::vector<Library> libraries;
using LibraryKey = std::pair<Library::Type, std::string>;
extern std::map<LibraryKey, std::size_t> library_dict;

//==============================================================================
// Non-member functions
//==============================================================================

extern "C" void read_cross_sections_xml();
void read_ce_cross_sections_xml();

} // namespace openmc

#endif // OPENMC_CROSS_SECTIONS_H
