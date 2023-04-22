#ifndef OPENMC_FILE_UTILS_H
#define OPENMC_FILE_UTILS_H

#include <fstream> // for ifstream
#include <string>
#include <sys/stat.h>

namespace openmc {

// TODO: replace with std::filesystem when switch to C++17 is made
//! Determine if a path is a directory
//! \param[in] path Path to check
//! \return Whether the path is a directory
inline bool dir_exists(const std::string& path)
{
  struct stat s;
  if (stat(path.c_str(), &s) != 0)
    return false;

  return s.st_mode & S_IFDIR;
}

//! Determine if a file exists
//! \param[in] filename Path to file
//! \return Whether file exists
inline bool file_exists(const std::string& filename)
{
  // rule out file being a path to a directory
  if (dir_exists(filename))
    return false;

  std::ifstream s {filename};
  return s.good();
}

} // namespace openmc

#endif // OPENMC_FILE_UTILS_H
