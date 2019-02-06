#include "openmc/error.h"

#include "openmc/message_passing.h"

#if defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h> // for isatty
#endif

#include <iomanip> // for setw
#include <iostream>

namespace openmc {

#ifdef OPENMC_MPI
void abort_mpi(int code)
{
  MPI_Abort(mpi::intracomm, code);
}
#endif

void warning(const std::string& message)
{
#ifdef _POSIX_VERSION
  // Make output yellow if user is in a terminal
  if (isatty(STDERR_FILENO)) {
    std::cerr << "\033[0;33m";
  }
#endif

  // Write warning at beginning
  std::cerr << " WARNING: ";

  // Set line wrapping and indentation
  int line_wrap = 80;
  int indent = 10;

  // Determine length of message
  int length = message.size();

  int i_start = 0;
  int line_len = line_wrap - indent + 1;
  while (i_start < length) {
    if (length - i_start < line_len) {
      // Remainder of message will fit on line
      std::cerr << message.substr(i_start) << '\n';
      break;

    } else {
      // Determine last space in current line
      std::string s = message.substr(i_start, line_len);
      auto pos = s.find_last_of(' ');

      // Write up to last space, or whole line if no space is present
      std::cerr << s.substr(0, pos) << '\n' << std::setw(indent) << " ";

      // Advance starting position
      i_start += (pos == std::string::npos) ? line_len : pos + 1;
    }
  }

#ifdef _POSIX_VERSION
  // Reset color for terminal
  if (isatty(STDERR_FILENO)) {
    std::cerr << "\033[0m";
  }
#endif
}

} // namespace openmc
