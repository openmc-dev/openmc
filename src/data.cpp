#include "openmc/data.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

std::string remove_space(std::string s)
{
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
  return s;
}

ZAI_data get_atomic_data(std::string line)
{
  ZAI_data atom_data;
  atom_data.mass_excess =
    std::stod(line.substr(30, 35) + "." + line.substr(36, 36));
  atom_data.binding_energy =
    std::stod(line.substr(56, 61) + "." + line.substr(62, 67));
  atom_data.mass =
    std::stod(line.substr(106, 109)) +
    1e-6 * std::stod(line.substr(110, 116) + '.' + line.substr(117, 123));
}

namespace openmc {

AtomicData::AtomicData(std::string data_file)
{
  // Open Data file
  std::ifstream myfile(data_file);
  // Skip header
  for (i = 0; i < line_in_header; i++) {
    std::string line;
    std::getline(myfile, line);
    std::string name = remove_space(line.substr(20, 22) + line.substr(16, 19));
    atomic_mass_data.insert(std::make_pair(name, get_atomic_data(line)));
  }
}

} // namespace openmc
