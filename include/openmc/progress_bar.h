#ifndef OPENMC_PROGRESSBAR_H
#define OPENMC_PROGRESSBAR_H

#include <string>

class ProgressBar {

public:  
  // Constructor
  ProgressBar();

  void set_value(double val);
  
private:
  std::string bar;
  char bar_old[72] = "???% |                                                                |";
  
};


#endif // OPENMC_PROGRESSBAR_H
  
