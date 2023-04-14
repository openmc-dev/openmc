#ifndef OPENMC_FILE_UTILS_H
#define OPENMC_FILE_UTILS_H

#include <algorithm> // any_of
#include <cctype>    // for isalpha
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

// Gets the file extension of whatever string is passed in. This is defined as
// a sequence of strictly alphanumeric characters which follow the last period,
// i.e. at least one alphabet character is present, and zero or more numbers.
// If such a sequence of characters is not found, an empty string is returned.
inline std::string get_file_extension(const std::string& filename)
{
  // check that at least one letter is present
  const std::string::size_type last_period_pos = filename.find_last_of('.');

  // no file extension
  if (last_period_pos == std::string::npos)
    return "";

  const std::string ending = filename.substr(last_period_pos + 1);

  // check that at least one character is present.
  const bool has_alpha = std::any_of(ending.begin(), ending.end(),
    [](char x) { return static_cast<bool>(std::isalpha(x)); });
  return has_alpha ? ending : "";
}

} // namespace openmc

#endif // OPENMC_FILE_UTILS_H
