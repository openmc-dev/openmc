#include "openmc/file_utils.h"

#include <filesystem>

namespace openmc {

bool dir_exists(const std::string& path)
{
  std::filesystem::path d(path);
  return std::filesystem::is_directory(d);
}

bool file_exists(const std::string& filename)
{
  std::filesystem::path p(filename);
  if (!std::filesystem::exists(p)) {
    return false;
  }
  if (std::filesystem::is_directory(p)) {
    return false;
  }
  return true;
}

std::string dir_name(const std::string& filename)
{
  std::filesystem::path p(filename);
  return (p.parent_path()).string();
}

std::string get_file_extension(const std::string& filename)
{
  std::filesystem::path p(filename);
  auto ext = p.extension();
  if (!ext.empty()) {
    // path::extension includes the period
    return ext.string().substr(1);
  }
  return "";
}

} // namespace openmc
