#ifndef ERROR_H
#define ERROR_H

#include <cstring>
#include <string>
#include <sstream>


namespace openmc {


extern "C" void fatal_error_from_c(const char* message, int message_len);
extern "C" void warning_from_c(const char* message, int message_len);


inline
void fatal_error(const char *message)
{
  fatal_error_from_c(message, strlen(message));
}

inline
void fatal_error(const std::string &message)
{
  fatal_error_from_c(message.c_str(), message.length());
}

inline
void fatal_error(const std::stringstream &message)
{
  fatal_error(message.str());
}

inline
void warning(const std::string& message)
{
  warning_from_c(message.c_str(), message.length());
}

inline
void warning(const std::stringstream& message)
{
  warning(message.str());
}

} // namespace openmc
#endif // ERROR_H
