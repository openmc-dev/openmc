#ifndef VARIANCE_REDUCTION_H
#define VARIANCE_REDUCTION_H 1

#include "pugixml/pugixml.hpp"
#include "hdf5.h"

namespace openmc {

  void read_importances(pugi::xml_node *node);

  void read_weight_windows(pugi::xml_node *node);

} // end openmc namespace

#endif
