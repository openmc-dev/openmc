#include "openmc/atomic_mass_data.h"

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

// Read a line from atomic mass data and pull usesull content
atomic_data read_atomic_data(std::string line)
{
  atomic_data atom_data;
  atom_data.mass_excess =
    std::stod(line.substr(29, 6) + "." + line.substr(36, 6));
  atom_data.binding_energy =
    std::stod(line.substr(55, 6) + "." + line.substr(62, 6));
  atom_data.mass =
    std::stod(line.substr(106, 3)) +
    1e-6 * std::stod(line.substr(110, 6) + '.' + line.substr(117, 3));
  return atom_data;
}

namespace openmc {
namespace data {

AtomicData::AtomicData(std::string data_file)
{
  // Open Data file
  std::ifstream myfile(data_file);
  if (!myfile) {
    std::cout << "NO file found in " << data_file << std::endl;
  }
  // Skip header
  for (int i = 0; i < line_in_header; i++) {
    std::string line;
    std::getline(myfile, line);
  }
  std::string line;
  std::getline(myfile, line);
  // read the data
  do {
    std::string name = remove_space(line.substr(20, 2) + line.substr(16, 3));
    atomic_data read = read_atomic_data(line);
    atomic_mass_data.insert(std::make_pair(name, read_atomic_data(line)));
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
