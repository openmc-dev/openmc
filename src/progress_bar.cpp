
#include "openmc/progress_bar.h"

#include <sstream>
#include <iostream>
#include <iomanip>

#ifdef UNIX
#include <unistd.h>
#endif

#define BAR_WIDTH 72

#ifdef UNIX
int check_isatty(int fd) {
  return isatty(fd);
}
#endif

bool is_terminal() {
#ifdef UNIX
  return check_isatty(STDOUT_FILENO) != 0;
#else
  return false;
#endif
}

ProgressBar::ProgressBar() {
  bar = "";
  set_value(0.0);
}

void
ProgressBar::set_value(double val) {

  if (!is_terminal()) return;
  
  // set the bar percentage
  if (val >= 100.0) {
    bar.append("100");
  } else if (val <= 0.0) {
    bar.append("  0");
  } else {
    std::stringstream ss;
    ss << std::setfill(' ') << std::setw(3) << (int)val;
    bar.append(ss.str());
  }

  bar.append("% |");

  // remaining width of the bar
  int remain = BAR_WIDTH - bar.size() - 2;
  
  // set the bar width
  if (val >= 100.0) {
    bar.append(remain, '=');
  } else {
    int width = (int)((double)remain*val/100);
    bar.append(width, '=');
    bar.append(1, '>');
    bar.append(remain-width-1, ' ');
  }

  bar.append("|+");
  
  // write the bar
  std::cout << '\r' << bar << std::flush;

  // reset the bar value
  bar = "";
}
