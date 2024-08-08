#ifndef OPENMC_STRING_UTILS_H
#define OPENMC_STRING_UTILS_H

#include <string>

#include "xtensor/xarray.hpp"

#include "openmc/vector.h"

namespace openmc {

std::string& strtrim(std::string& s);

char* strtrim(char* c_str);

std::string to_element(const std::string& name);

void to_lower(std::string& str);

int word_count(const std::string& str);

vector<std::string> split(const std::string& in);

bool ends_with(const std::string& value, const std::string& ending);

bool starts_with(const std::string& value, const std::string& beginning);

template<typename T>
inline std::string concatenate(const T& values, const std::string& del = " ")
{
  std::ostringstream oss;
  auto it = values.begin();
  if (it != values.end()) {
    oss << *it; // Insert the first element without a delimiter
    ++it;
  }
  for (; it != values.end(); ++it) {
    oss << del << *it;
  }
  return oss.str();
}

} // namespace openmc
#endif // OPENMC_STRING_UTILS_H
