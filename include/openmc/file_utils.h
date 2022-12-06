#ifndef OPENMC_FILE_UTILS_H
#define OPENMC_FILE_UTILS_H

#include <fstream> // for ifstream
#include <string>
#include <sys/stat.h>

namespace openmc {

// TODO: replace with std::filesysem when switch to C++17 is made
//! Determine if a path is a directory
//! \param[in] path Path to check
//! \return Whether the path is a directory
inline bool is_dir(const std::string& path) {
  struct stat s;
  if (stat(path.c_str(), &s) != 0) return false;

  return s.st_mode & S_IFDIR;
}

//! Determine if a file exists
//! \param[in] filename Path to file
//! \return Whether file exists
inline bool file_exists(const std::string& filename)
{
  // rule out file being a directory path
  if (is_dir(filename)) return false;

  std::ifstream s {filename};
  s.seekg(0, std::ios::beg);
  return s.good();
}

} // namespace openmc

#endif // OPENMC_FILE_UTILS_H
