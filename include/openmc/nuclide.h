#ifndef OPENMC_NUCLIDE_H
#define OPENMC_NUCLIDE_H

#include "openmc/constants.h"

namespace openmc {

//===============================================================================
//! Cached microscopic cross sections for a particular nuclide at the current
//! energy
//===============================================================================

struct NuclideMicroXS {
  // Microscopic cross sections in barns
  double total;            //!< total cross section
  double absorption;       //!< absorption (disappearance)
  double fission;          //!< fission
  double nu_fission;       //!< neutron production from fission

  double elastic;          //!< If sab_frac is not 1 or 0, then this value is
                           //!<   averaged over bound and non-bound nuclei
  double thermal;          //!< Bound thermal elastic & inelastic scattering
  double thermal_elastic;  //!< Bound thermal elastic scattering
  double photon_prod;      //!< microscopic photon production xs

  // Cross sections for depletion reactions (note that these are not stored in
  // macroscopic cache)
  double reaction[DEPLETION_RX.size()];

  // Indicies and factors needed to compute cross sections from the data tables
  int index_grid;        //!< Index on nuclide energy grid
  int index_temp;        //!< Temperature index for nuclide
  double interp_factor;  //!< Interpolation factor on nuc. energy grid
  int index_sab {-1};    //!< Index in sab_tables
  int index_temp_sab;    //!< Temperature index for sab_tables
  double sab_frac;       //!< Fraction of atoms affected by S(a,b)
  bool use_ptable;       //!< In URR range with probability tables?

  // Energy and temperature last used to evaluate these cross sections.  If
  // these values have changed, then the cross sections must be re-evaluated.
  double last_E {0.0};      //!< Last evaluated energy
  double last_sqrtkT {0.0}; //!< Last temperature in sqrt(Boltzmann constant
                            //!< * temperature (eV))
};

//===============================================================================
// MATERIALMACROXS contains cached macroscopic cross sections for the material a
// particle is traveling through
//===============================================================================

  struct MaterialMacroXS {
    double total;         //!< macroscopic total xs
    double absorption;    //!< macroscopic absorption xs
    double fission;       //!< macroscopic fission xs
    double nu_fission;    //!< macroscopic production xs
    double photon_prod;   //!< macroscopic photon production xs

    // Photon cross sections
    double coherent;        //!< macroscopic coherent xs
    double incoherent;      //!< macroscopic incoherent xs
    double photoelectric;   //!< macroscopic photoelectric xs
    double pair_production; //!< macroscopic pair production xs
  };

} // namespace openmc

#endif // OPENMC_NUCLIDE_H
