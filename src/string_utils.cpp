#include "openmc/string_utils.h"

#include <algorithm> // for equal
#include <cctype>    // for tolower, isspace
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


std::string to_element(const std::string& name) {
  int pos = name.find_first_of("0123456789");
  return name.substr(0, pos);
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
  while (stream >> dum) {count++;}
  return count;
}

std::vector<std::string> split(const std::string& in)
{
  std::vector<std::string> out;

  for (int i = 0; i < in.size(); ) {
    // Increment i until we find a non-whitespace character.
    if (std::isspace(in[i])) {
      i++;

    } else {
      // Find the next whitespace character at j.
      int j = i + 1;
      while (j < in.size() && std::isspace(in[j]) == 0) {j++;}

      // Push-back everything between i and j.
      out.push_back(in.substr(i, j-i));
      i = j + 1; // j is whitespace so leapfrog to j+1
    }
  }

  return out;
}

bool ends_with(const std::string& value, const std::string& ending)
{
  if (ending.size() > value.size()) return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

bool starts_with(const std::string& value, const std::string& beginning)
{
  if (beginning.size() > value.size()) return false;
  return std::equal(beginning.begin(), beginning.end(), value.begin());
}

} // namespace openmc
