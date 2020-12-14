#include "openmc/summary.h"

#include "openmc/cell.h"
#include "openmc/hdf5_interface.h"
#include "openmc/lattice.h"
#include "openmc/material.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/output.h"
#include "openmc/surface.h"
#include "openmc/settings.h"

namespace openmc {

void write_summary()
{
  // Display output message
  write_message("Writing summary.h5 file...", 5);

  // Create a new file using default properties.
  hid_t file = file_open("summary.h5", 'w');

  write_header(file);
  write_nuclides(file);
  write_geometry(file);
  write_materials(file);

  // Terminate access to the file.
  file_close(file);
}

void write_header(hid_t file)
{
  // Write filetype and version info
  write_attribute(file, "filetype", "summary");
  write_attribute(file, "version", VERSION_SUMMARY);
  write_attribute(file, "openmc_version", VERSION);
#ifdef GIT_SHA1
  write_attribute(file, "git_sha1", GIT_SHA1);
#endif

  // Write current date and time
  write_attribute(file, "date_and_time", time_stamp());
}

void write_nuclides(hid_t file)
{
  // Build vectors of nuclide names and awrs while only sorting nuclides from
  // macroscopics
  std::vector<std::string> nuc_names;
  std::vector<std::string> macro_names;
  std::vector<double> awrs;

  for (int i = 0; i < data::nuclides.size(); ++i) {
    if (settings::run_CE) {
      const auto& nuc {data::nuclides[i]};
      nuc_names.push_back(nuc->name_);
      awrs.push_back(nuc->awr_);
    } else {
      const auto& nuc {data::mg.nuclides_[i]};
      if (nuc.awr != MACROSCOPIC_AWR) {
        nuc_names.push_back(nuc.name);
        awrs.push_back(nuc.awr);
      } else {
        macro_names.push_back(nuc.name);
      }
    }
  }

  hid_t nuclide_group = create_group(file, "nuclides");
  write_attribute(nuclide_group, "n_nuclides", nuc_names.size());
  hid_t macro_group = create_group(file, "macroscopics");
  write_attribute(macro_group, "n_macroscopics", macro_names.size());
  // Write nuclide names and awrs
  if (!nuc_names.empty()) {
    // Write useful data from nuclide objects
    write_dataset(nuclide_group, "names", nuc_names);
    write_dataset(nuclide_group, "awrs", awrs);
  }
  if (!macro_names.empty()) {
    // Write useful data from macroscopic objects
    write_dataset(macro_group, "names", macro_names);
  }
  close_group(nuclide_group);
  close_group(macro_group);
}

void write_geometry(hid_t file)
{
  auto geom_group = create_group(file, "geometry");

#ifdef DAGMC
  if (settings::dagmc) {
    write_attribute(geom_group, "dagmc", 1);
    close_group(geom_group);
    return;
  }
#endif

  write_attribute(geom_group, "n_cells", model::cells.size());
  write_attribute(geom_group, "n_surfaces", model::surfaces.size());
  write_attribute(geom_group, "n_universes", model::universes.size());
  write_attribute(geom_group, "n_lattices", model::lattices.size());

  auto cells_group = create_group(geom_group, "cells");
  for (const auto& c : model::cells) c->to_hdf5(cells_group);
  close_group(cells_group);

  auto surfaces_group = create_group(geom_group, "surfaces");
  for (const auto& surf : model::surfaces) surf->to_hdf5(surfaces_group);
  close_group(surfaces_group);

  auto universes_group = create_group(geom_group, "universes");
  for (const auto& u : model::universes) u->to_hdf5(universes_group);
  close_group(universes_group);

  auto lattices_group = create_group(geom_group, "lattices");
  for (const auto& lat : model::lattices) lat->to_hdf5(lattices_group);
  close_group(lattices_group);

  close_group(geom_group);
}

void write_materials(hid_t file)
{
  // write number of materials
  write_dataset(file, "n_materials", model::materials.size());

  hid_t materials_group = create_group(file, "materials");
  for (const auto& mat : model::materials) {
    mat->to_hdf5(materials_group);
  }
  close_group(materials_group);
}

} // namespace openmc
