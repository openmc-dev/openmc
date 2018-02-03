#ifndef ERROR_H
#define ERROR_H

#include <cstring>
#include <string>
#include <sstream>


namespace openmc {


extern "C" void fatal_error_from_c(const char *message, int message_len);


void fatal_error(const char *message)
{
  fatal_error_from_c(message, strlen(message));
}


void fatal_error(const std::string &message)
{
  fatal_error_from_c(message.c_str(), message.length());
}


void fatal_error(const std::stringstream &message)
{
  std::string out {message.str()};
  fatal_error_from_c(out.c_str(), out.length());
}

} // namespace openmc
#endif // ERROR_H
