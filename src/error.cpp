#include "openmc/error.h"

#include "openmc/message_passing.h"
#include "openmc/settings.h"

#if defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h> // for isatty
#endif

#include <cstdlib> // for exit
#include <iomanip> // for setw
#include <iostream>

//==============================================================================
// Global variables / constants
//==============================================================================

// Error codes
int OPENMC_E_UNASSIGNED {-1};
int OPENMC_E_ALLOCATE {-2};
int OPENMC_E_OUT_OF_BOUNDS {-3};
int OPENMC_E_INVALID_SIZE {-4};
int OPENMC_E_INVALID_ARGUMENT {-5};
int OPENMC_E_INVALID_TYPE {-6};
int OPENMC_E_INVALID_ID {-7};
int OPENMC_E_GEOMETRY {-8};
int OPENMC_E_DATA {-9};
int OPENMC_E_PHYSICS {-10};
int OPENMC_E_WARNING {1};

// Error message
char openmc_err_msg[256];

//==============================================================================
// Functions
//==============================================================================

namespace openmc {

#ifdef OPENMC_MPI
void abort_mpi(int code)
{
  MPI_Abort(mpi::intracomm, code);
}
#endif

void output(const std::string& message, std::ostream& out, int indent)
{
  // Set line wrapping and indentation
  int line_wrap = 80;

  // Determine length of message
  int length = message.size();

  int i_start = 0;
  int line_len = line_wrap - indent + 1;
  while (i_start < length) {
    if (length - i_start < line_len) {
      // Remainder of message will fit on line
      out << message.substr(i_start) << std::endl;
      break;

    } else {
      // Determine last space in current line
      std::string s = message.substr(i_start, line_len);
      auto pos = s.find_last_of(' ');

      // Write up to last space, or whole line if no space is present
      out << s.substr(0, pos) << '\n' << std::setw(indent) << " ";

      // Advance starting position
      i_start += (pos == std::string::npos) ? line_len : pos + 1;
    }
  }
}

void warning(const std::string& message)
{
#ifdef _POSIX_VERSION
  // Make output yellow if user is in a terminal
  if (isatty(STDERR_FILENO)) {
    std::cerr << "\033[0;33m";
  }
#endif

  // Write warning
  std::cerr << " WARNING: ";
  output(message, std::cerr, 10);

#ifdef _POSIX_VERSION
  // Reset color for terminal
  if (isatty(STDERR_FILENO)) {
    std::cerr << "\033[0m";
  }
#endif
}

void write_message(const std::string& message, int level)
{
  // Only allow master to print to screen
  if (!mpi::master) return;

  if (level <= settings::verbosity) {
    std::cout << " ";
    output(message, std::cout, 1);
  }
}

void fatal_error(const std::string& message, int err)
{
#ifdef _POSIX_VERSION
  // Make output red if user is in a terminal
  if (isatty(STDERR_FILENO)) {
    std::cerr << "\033[0;31m";
  }
#endif

  // Write error message
  std::cerr << " ERROR: ";
  output(message, std::cerr, 8);

#ifdef _POSIX_VERSION
  // Reset color for terminal
  if (isatty(STDERR_FILENO)) {
    std::cerr << "\033[0m";
  }
#endif

#ifdef OPENMC_MPI
  MPI_Abort(mpi::intracomm, err);
#endif

  // Abort the program
  std::exit(err);
}

} // namespace openmc
