#ifndef OPENMC_ERROR_H
#define OPENMC_ERROR_H

#include <cstring>
#include <string>
#include <sstream>

#include "openmc/capi.h"

namespace openmc {


extern "C" void fatal_error_from_c(const char* message, int message_len);
extern "C" void warning_from_c(const char* message, int message_len);
extern "C" void write_message_from_c(const char* message, int message_len,
                                     int level);

inline void
set_errmsg(const char* message)
{
  std::strcpy(openmc_err_msg, message);
}

inline void
set_errmsg(const std::string& message)
{
  std::strcpy(openmc_err_msg, message.c_str());
}

inline void
set_errmsg(const std::stringstream& message)
{
  std::strcpy(openmc_err_msg, message.str().c_str());
}

inline
void fatal_error(const char* message)
{
  fatal_error_from_c(message, std::strlen(message));
}

inline
void fatal_error(const std::string& message)
{
  fatal_error_from_c(message.c_str(), message.length());
}

inline
void fatal_error(const std::stringstream& message)
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

inline
void write_message(const char* message, int level)
{
  write_message_from_c(message, std::strlen(message), level);
}

inline
void write_message(const std::string& message, int level)
{
  write_message_from_c(message.c_str(), message.length(), level);
}

inline
void write_message(const std::stringstream& message, int level)
{
  write_message(message.str(), level);
}

} // namespace openmc
#endif // OPENMC_ERROR_H
