#ifndef OPENMC_FILE_UTILS_H
#define OPENMC_FILE_UTILS_H

#include <fstream> // for ifstream
#include <string>

namespace openmc {

inline bool file_exists(const std::string& filename)
{
  std::ifstream s {filename};
  return s.good();
}

}

#endif // OPENMC_FILE_UTILS_H