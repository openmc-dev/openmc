#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <string>
#include <vector>


namespace openmc {

std::vector<std::string>
split(const std::string &in)
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

} // namespace openmc
#endif // STRING_UTILS_H
