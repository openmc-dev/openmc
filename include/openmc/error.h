#ifndef OPENMC_ERROR_H
#define OPENMC_ERROR_H

#include <cstring>
#include <sstream>
#include <string>

#include <fmt/format.h>

#include "openmc/capi.h"
#include "openmc/settings.h"

#if defined(__GNUC__) || defined(__clang__)
#define UNREACHABLE() __builtin_unreachable()
#else
#define UNREACHABLE() (void)0
#endif

namespace openmc {

void set_errmsg(const char* message);
void set_errmsg(const std::string& message);
void set_errmsg(const std::stringstream& message);
const char* get_errmsg();

[[noreturn]] void fatal_error(const std::string& message, int err = -1);

[[noreturn]] inline void fatal_error(const std::stringstream& message)
{
  fatal_error(message.str());
}

[[noreturn]] inline void fatal_error(const char* message)
{
  fatal_error(std::string {message, std::strlen(message)});
}

void warning(const std::string& message);

inline void warning(const std::stringstream& message)
{
  warning(message.str());
}

void write_message(const std::string& message, int level = 0);

inline void write_message(const std::stringstream& message, int level)
{
  write_message(message.str(), level);
}

template<typename... Params>
void write_message(
  int level, const std::string& message, const Params&... fmt_args)
{
  if (settings::verbosity >= level) {
    write_message(fmt::format(fmt::runtime(message), fmt_args...));
  }
}

template<typename... Params>
void write_message(const std::string& message, const Params&... fmt_args)
{
  write_message(fmt::format(fmt::runtime(message), fmt_args...));
}

#ifdef OPENMC_MPI
extern "C" void abort_mpi(int code);
#endif

} // namespace openmc
#endif // OPENMC_ERROR_H
