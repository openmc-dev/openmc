#include "openmc/summary.h"

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/lattice.h"
#include "openmc/material.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/output.h"
#include "openmc/settings.h"
#include "openmc/surface.h"

namespace openmc {

void write_summary()
{
  // Display output message
  write_message("Writing summary.h5 file...", 5);

  // Set filename for summary file
  std::string filename = fmt::format("{}summary.h5", settings::path_output);

  // Create a new file using default properties.
  hid_t file = file_open(filename, 'w');

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
  vector<std::string> nuc_names;
  vector<std::string> macro_names;
  vector<double> awrs;

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

  write_attribute(geom_group, "n_cells", model::cells.size());
  write_attribute(geom_group, "n_surfaces", model::surfaces.size());
  write_attribute(geom_group, "n_universes", model::universes.size());
  write_attribute(geom_group, "n_lattices", model::lattices.size());

  auto cells_group = create_group(geom_group, "cells");
  for (const auto& c : model::cells)
    c->to_hdf5(cells_group);
  close_group(cells_group);

  auto surfaces_group = create_group(geom_group, "surfaces");
  for (const auto& surf : model::surfaces)
    surf->to_hdf5(surfaces_group);
  close_group(surfaces_group);

  auto universes_group = create_group(geom_group, "universes");
  for (const auto& u : model::universes)
    u->to_hdf5(universes_group);
  close_group(universes_group);

  auto lattices_group = create_group(geom_group, "lattices");
  for (const auto& lat : model::lattices)
    lat->to_hdf5(lattices_group);
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

//==============================================================================
// C API
//==============================================================================

extern "C" int openmc_properties_export(const char* filename)
{
  // Only write from master process
  if (!mpi::master)
    return 0;

  // Set a default filename if none was passed
  std::string name = filename ? filename : "properties.h5";

  // Display output message
  auto msg = fmt::format("Exporting properties to {}...", name);
  write_message(msg, 5);

  // Create a new file using default properties.
  hid_t file = file_open(name, 'w');

  // Write metadata
  write_attribute(file, "filetype", "properties");
  write_attribute(file, "version", VERSION_STATEPOINT);
  write_attribute(file, "openmc_version", VERSION);
#ifdef GIT_SHA1
  write_attribute(file, "git_sha1", GIT_SHA1);
#endif
  write_attribute(file, "date_and_time", time_stamp());
  write_attribute(file, "path", settings::path_input);

  // Write cell properties
  auto geom_group = create_group(file, "geometry");
  write_attribute(geom_group, "n_cells", model::cells.size());
  auto cells_group = create_group(geom_group, "cells");
  for (const auto& c : model::cells) {
    c->export_properties_hdf5(cells_group);
  }
  close_group(cells_group);
  close_group(geom_group);

  // Write material properties
  hid_t materials_group = create_group(file, "materials");
  write_attribute(materials_group, "n_materials", model::materials.size());
  for (const auto& mat : model::materials) {
    mat->export_properties_hdf5(materials_group);
  }
  close_group(materials_group);

  // Terminate access to the file.
  file_close(file);
  return 0;
}

extern "C" int openmc_properties_import(const char* filename)
{
  // Display output message
  auto msg = fmt::format("Importing properties from {}...", filename);
  write_message(msg, 5);

  // Create a new file using default properties.
  if (!file_exists(filename)) {
    set_errmsg(fmt::format("File '{}' does not exist.", filename));
    return OPENMC_E_INVALID_ARGUMENT;
  }
  hid_t file = file_open(filename, 'r');

  // Ensure the filetype is correct
  std::string filetype;
  read_attribute(file, "filetype", filetype);
  if (filetype != "properties") {
    file_close(file);
    set_errmsg(fmt::format("File '{}' is not a properties file.", filename));
    return OPENMC_E_INVALID_ARGUMENT;
  }

  // Make sure number of cells matches
  auto geom_group = open_group(file, "geometry");
  int32_t n;
  read_attribute(geom_group, "n_cells", n);
  if (n != openmc::model::cells.size()) {
    close_group(geom_group);
    file_close(file);
    set_errmsg(fmt::format(
      "Number of cells in {} doesn't match current model.", filename));
    return OPENMC_E_GEOMETRY;
  }

  // Read cell properties
  auto cells_group = open_group(geom_group, "cells");
  try {
    for (const auto& c : model::cells) {
      c->import_properties_hdf5(cells_group);
    }
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_UNASSIGNED;
  }
  close_group(cells_group);
  close_group(geom_group);

  // Make sure number of cells matches
  auto materials_group = open_group(file, "materials");
  read_attribute(materials_group, "n_materials", n);
  if (n != openmc::model::materials.size()) {
    close_group(materials_group);
    file_close(file);
    set_errmsg(fmt::format(
      "Number of materials in {} doesn't match current model.", filename));
    return OPENMC_E_GEOMETRY;
  }

  // Read material properties
  for (const auto& mat : model::materials) {
    mat->import_properties_hdf5(materials_group);
  }
  close_group(materials_group);

  // Terminate access to the file.
  file_close(file);
  return 0;
}

} // namespace openmc
