#include "openmc/atomic_mass_data.h"
#include "openmc/error.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// Clean string to remove unrequired spaces
std::string remove_space(std::string s)
{
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
  return s;
}

namespace openmc {
namespace data {

// Read a line from atomic mass data and pull usesull content and return a
// pair<nuc_name, atomic_mass_data>
std::pair<std::string, atomic_data> read(std::string line)
{
  std::string nuc_name = remove_space(line.substr(20, 2) + line.substr(16, 3));
  atomic_data atom_data;
  atom_data.mass_excess =
    std::stod(line.substr(29, 6) + "." + line.substr(36, 6));
  atom_data.binding_energy =
    std::stod(line.substr(55, 6) + "." + line.substr(62, 6));
  atom_data.mass =
    std::stod(line.substr(106, 3)) +
    1e-6 * std::stod(line.substr(110, 6) + '.' + line.substr(117, 3));
  return std::make_pair(nuc_name, atom_data);
}

AtomicData::AtomicData(std::string data_file)
{
  // Open Data file
  std::ifstream myfile(data_file);
  if (!myfile) {
    fatal_error("Atomic mass data file '" + data_file + "' does not exist.");
  }
  std::string line = "";

  // Skip header and read the first line with data
  for (int i = 0; i < line_in_header + 1; i++) {
    std::getline(myfile, line);
  }

  do {
    using openmc::data::read;
    atomic_mass_data.insert(read(line));
    std::getline(myfile, line);
  } while (!myfile.eof());
}

double AtomicData::get_atomic_mass(std::string nuclide) const
{
  return get_atomic_data(nuclide).mass;
}

double AtomicData::get_atomic_mass_excess(std::string nuclide) const
{
  return get_atomic_data(nuclide).mass_excess;
}

double AtomicData::get_atomic_binding_energy(std::string nuclide) const
{
  return get_atomic_data(nuclide).binding_energy;
}

atomic_data AtomicData::get_atomic_data(std::string nuclide) const
{
  auto pos = atomic_mass_data.find(nuclide);

  if (pos == atomic_mass_data.end()) {
    throw std::out_of_range(
      "Nuclide " + nuclide + " not in the atomic mass data.");
  } else {
    return pos->second;
  }
}
} // namespace data
} // namespace openmc
