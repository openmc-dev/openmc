#ifndef OPENMC_FILE_UTILS_H
#define OPENMC_FILE_UTILS_H

#include <fstream> // for ifstream
#include <string>

namespace openmc {

//! Determine if a file exists
//! \param[in] filename Path to file
//! \return Whether file exists
inline bool file_exists(const std::string& filename)
{
  std::ifstream s {filename};
  return s.good();
}

} // namespace openmc

#endif // OPENMC_FILE_UTILS_H
