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
  Type type_; //!< Type of data library
  std::vector<std::string> materials_; //!< Materials contained in library
  std::string path_; //!< File path to library
};

using LibraryKey = std::pair<Library::Type, std::string>;

//==============================================================================
// Global variable declarations
//==============================================================================

namespace data {

//!< Data libraries
extern std::vector<Library> libraries;

//! Maps (type, name) to index in libraries
extern std::map<LibraryKey, std::size_t> library_map;

} // namespace data

//==============================================================================
// Non-member functions
//==============================================================================

//! Read cross sections file (either XML or multigroup H5) and populate data
//! libraries
void read_cross_sections_xml();


//! Load nuclide and thermal scattering data from HDF5 files
//
//! \param[in] nuc_temps Temperatures for each nuclide in [K]
//! \param[in] thermal_temps Temperatures for each thermal scattering table in [K]
void read_ce_cross_sections(const std::vector<std::vector<double>>& nuc_temps,
  const std::vector<std::vector<double>>& thermal_temps);

//! Read cross_sections.xml and populate data libraries
void read_ce_cross_sections_xml();

void library_clear();

} // namespace openmc

#endif // OPENMC_CROSS_SECTIONS_H
