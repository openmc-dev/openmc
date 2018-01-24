#ifndef ERROR_H
#define ERROR_H

#include <cstring>
#include <string>


extern "C" void fatal_error_from_c(const char *message, int message_len);


void fatal_error(const char *message)
{
  fatal_error_from_c(message, strlen(message));
}


void fatal_error(const std::string &message)
{
  fatal_error_from_c(message.c_str(), message.length());
}

#endif // ERROR_H
