#include "openmc/variance_reduction.h"

#include <iostream>
#include <array>
#include <cmath>
#include <sstream>

#include "error.h"

namespace openmc {
namespace variance_reduction {
bool importance_splitting = false;
  bool weight_splitting = false;
  std::vector<double> importances;
extern "C" void
read_variance_reduction(pugi::xml_node *node)
{
  int n_imp = 0; // number of importance blocks
  int n_ww = 0; // number of weight window blocks
  for (pugi::xml_node imp_node: node->children("importance")) {n_imp++;}
  for (pugi::xml_node ww_node: node->children("weight_window")) {n_ww++;}

  // vr not required
  if (n_imp == 0 && n_ww == 0 ) {
    std::cout << "No VR set - analogue mode" << std::endl;
  } else if ( n_imp >= 1 && n_ww >= 1 ) { // too many 
    std::stringstream err_msg;
    err_msg << "importances and weight windows are input "
               "one or ther other is allowed ";
    fatal_error(err_msg);
  }
  
  // read the importance node
  if (n_imp == 1) {
    read_importances(node);
  } else if ( n_ww == 1 ) {
    //read_weight_windows(node);
  } 
  
}

// read the importances
void
read_importances(pugi::xml_node *node)
{
  std::cout << "Found imp" << std::endl;
}

} // namespace variance_reduction 
} // namespace openmc
