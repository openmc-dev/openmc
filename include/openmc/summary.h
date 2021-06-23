#ifndef OPENMC_SUMMARY_H
#define OPENMC_SUMMARY_H

#include <hdf5.h>

namespace openmc {

void write_summary();
void write_header(hid_t file);
void write_nuclides(hid_t file);
void write_geometry(hid_t file);
void write_materials(hid_t file);

}

#endif // OPENMC_SUMMARY_H
