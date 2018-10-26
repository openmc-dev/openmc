#include "openmc/string_functions.h"
#include <sstream>

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

int word_count(std::string const& str)
{
  std::stringstream stream(str);
  std::string dum;
  int count = 0;
  while(stream >> dum) { count++; }
  return count;
}

} // namespace openmc
