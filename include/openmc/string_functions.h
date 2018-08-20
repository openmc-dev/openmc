//! \file string_functions.h
//! A collection of helper routines for C-strings and STL strings

#ifndef OPENMC_STRING_FUNCTIONS_H
#define OPENMC_STRING_FUNCTIONS_H

#include <string>

namespace openmc {

std::string& strtrim(std::string& s);

char* strtrim(char* c_str);

void to_lower(std::string& str);

} // namespace openmc
#endif // STRING_FUNCTIONS_H
