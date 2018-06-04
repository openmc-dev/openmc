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

void strtrim(char* str)
{
    int start = 0; // number of leading spaces
    char* buffer = str;

    while (*str && *str++ == ' ') ++start;

    while (*str++); // move to end of string

    // backup over trailing spaces
    int end = str - buffer - 1;
    while (end > 0 && buffer[end - 1] == ' ') --end;
    buffer[end] = 0; // remove trailing spaces

    // exit if no leading spaces or string is now empty
    if (end <= start || start == 0) return;
    str = buffer + start;

    while ((*buffer++ = *str++));  // remove leading spaces: K&R
}

std::string strtrim(std::string in_str)
{
  std::string str = in_str;
  // perform the left trim
  str.erase(str.begin(), std::find_if(str.begin(), str.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
  // perform the right trim
  str.erase(std::find_if(str.rbegin(), str.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
                      str.end());
}

void to_lower(std::string& str)
{
  for (int i = 0; i < str.size(); i++) str[i] = std::tolower(str[i]);
}

} // namespace openmc
#endif // STRING_FUNCTIONS_H