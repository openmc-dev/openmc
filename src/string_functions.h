//! \file string_functions.h
//! A collection of helper routines for C-strings and STL strings

#ifndef STRING_FUNCTIONS_H
#define STRING_FUNCTIONS_H

// for string functions
#include <functional>
#include <cctype>
#include <locale>
#include <string>

namespace openmc {

std::string& strtrim(std::string& s)
{
  const char* t = " \t\n\r\f\v";
  s.erase(s.find_last_not_of(t) + 1);
  s.erase(0, s.find_first_not_of(t));
  return s;
}


char* strtrim(char* c_str)
{
  std::string std_str;
  std_str.assign(c_str);
  strtrim(std_str);
  int length = std_str.copy(c_str, std_str.size());
  c_str[length] = '\0';
  return c_str;
}


void to_lower(std::string& str)
{
  for (int i = 0; i < str.size(); i++) str[i] = std::tolower(str[i]);
}

} // namespace openmc
#endif // STRING_FUNCTIONS_H