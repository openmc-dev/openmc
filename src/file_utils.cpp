#include "openmc/file_utils.h"

#include <algorithm> // any_of
#include <cctype>    // for isalpha
#include <fstream>   // for ifstream
#include <sys/stat.h>

namespace openmc {

bool dir_exists(const std::string& path)
{
  struct stat s;
  if (stat(path.c_str(), &s) != 0)
    return false;

  return s.st_mode & S_IFDIR;
}

bool file_exists(const std::string& filename)
{
  // rule out file being a path to a directory
  if (dir_exists(filename))
    return false;

  std::ifstream s {filename};
  return s.good();
}

std::string get_file_extension(const std::string& filename)
{
  // check that at least one letter is present
  const auto last_period_pos = filename.find_last_of('.');

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
