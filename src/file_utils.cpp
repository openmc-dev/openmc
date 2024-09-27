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

std::string dir_name(const std::string& filename)
{
  size_t pos = filename.find_last_of("\\/");
  return (std::string::npos == pos) ? "" : filename.substr(0, pos + 1);
}

std::string get_file_extension(const std::string& filename)
{
  // try our best to work on windows...
#if defined(_WIN32) || defined(_WIN64)
  const char sep_char = '\\';
#else
  const char sep_char = '/';
#endif

  // check that at least one letter is present
  const auto last_period_pos = filename.find_last_of('.');
  const auto last_sep_pos = filename.find_last_of(sep_char);

  // no file extension. In the first case, we are only given
  // a file name. In the second, we have been given a file path.
  // If that's the case, periods are allowed in directory names,
  // but have the interpretation as preceding a file extension
  // after the last separator.
  if (last_period_pos == std::string::npos ||
      (last_sep_pos < std::string::npos && last_period_pos < last_sep_pos))
    return "";

  const std::string ending = filename.substr(last_period_pos + 1);

  // check that at least one character is present.
  const bool has_alpha = std::any_of(ending.begin(), ending.end(),
    [](char x) { return static_cast<bool>(std::isalpha(x)); });
  return has_alpha ? ending : "";
}

} // namespace openmc
