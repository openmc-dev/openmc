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

std::string concatenate_vec(const std::vector<double>& temps);

std::string concatenate_xt(const xt::xarray<double>& temps);

} // namespace openmc
#endif // OPENMC_STRING_UTILS_H
