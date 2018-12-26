#include "openmc/wmp.h"

#include "openmc/hdf5_interface.h"

namespace openmc {

WindowedMultipole::WindowedMultipole(hid_t group)
{
  // Get name of nuclide from group, removing leading '/'
  name_ = object_name(group).substr(1);

  // Read scalar values.
  read_dataset(group, "spacing", spacing_);
  read_dataset(group, "sqrtAWR", sqrt_awr_);
  read_dataset(group, "E_min", E_min_);
  read_dataset(group, "E_max", E_max_);

  // Read the "data" array.  Use its shape to figure out the number of poles
  // and residue types in this data.
  read_dataset(group, "data", data_);
  int n_residues = data_.shape()[1] - 1;

  // Check to see if this data includes fission residues.
  fissionable_ = (n_residues == 3);

  // Read the "windows" array and use its shape to figure out the number of
  // windows.
  read_dataset(group, "windows", windows_);
  int n_windows = windows_.shape()[0];

  // Read the "broaden_poly" arrays.
  read_dataset(group, "broaden_poly", broaden_poly_);
  if (n_windows != broaden_poly_.shape()[0]) {
    fatal_error("broaden_poly array shape is not consistent with the windows "
      "array shape in WMP library for " + name_ + ".");
  }

  // Read the "curvefit" array.
  read_dataset(group, "curvefit", curvefit_);
  if (n_windows != broaden_poly_.shape()[0]) {
    fatal_error("curvefit array shape is not consistent with the windows "
      "array shape in WMP library for " + name_ + ".");
  }
  fit_order_ = curvefit_.shape()[1] - 1;
}

} // namespace openmc
