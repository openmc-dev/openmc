#ifndef OPENMC_CROSS_SECTIONS_H
#define OPENMC_CROSS_SECTIONS_H

#include "pugixml.hpp"

#include <map>
#include <string>

#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Library class
//==============================================================================

class Library {
public:
  // Types, enums
  enum class Type {
    neutron = 1,
    photon = 3,
    thermal = 2,
    multigroup = 4,
    wmp = 5
  };

  // Constructors
  Library() {};
  Library(pugi::xml_node node, const std::string& directory);

  // Comparison operator (for using in map)
  bool operator<(const Library& other) { return path_ < other.path_; }

  // Data members
  Type type_;                     //!< Type of data library
  vector<std::string> materials_; //!< Materials contained in library
  std::string path_;              //!< File path to library
};

using LibraryKey = std::pair<Library::Type, std::string>;

//==============================================================================
// Global variable declarations
//==============================================================================

namespace data {

//! Maps (type, name) to index in libraries
extern std::map<LibraryKey, std::size_t> library_map;

//!< Data libraries
extern vector<Library> libraries;

} // namespace data

//==============================================================================
// Non-member functions
//==============================================================================

//! Read cross sections file (either XML or multigroup H5) and populate data
//! libraries
void read_cross_sections_xml();

//! Read cross sections file (either XML or multigroup H5) and populate data
//! libraries
//! \param[in] root node of the cross_sections.xml
void read_cross_sections_xml(pugi::xml_node root);

//! Load nuclide and thermal scattering data from HDF5 files
//
//! \param[in] nuc_temps Temperatures for each nuclide in [K]
//! \param[in] thermal_temps Temperatures for each thermal scattering table in
//! [K]
void read_ce_cross_sections(const vector<vector<double>>& nuc_temps,
  const vector<vector<double>>& thermal_temps);

//! Read cross_sections.xml and populate data libraries
void read_ce_cross_sections_xml();

//! Load nuclide and thermal scattering data
void finalize_cross_sections();

void library_clear();

} // namespace openmc

#endif // OPENMC_CROSS_SECTIONS_H
