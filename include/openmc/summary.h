#ifndef OPENMC_SUMMARY_H
#define OPENMC_SUMMARY_H

#include <hdf5.h>

namespace openmc {

void write_summary();
void write_header(hid_t file);
void write_nuclides(hid_t file);
void write_geometry(hid_t file);
void write_materials(hid_t file);

//! Export physical properties for model
//! \param[in] filename Filename to write to
//! \return Error code
extern "C" int openmc_properties_export(const char* filename);

//! Import physical properties for model
//! \param[in] filename Filename to read from
// \return Error code
extern "C" int openmc_properties_import(const char* filename);

}

#endif // OPENMC_SUMMARY_H
