#include "openmc/tallies/filter.h"

#include <string>

#include "openmc/capi.h"
#include "openmc/constants.h"  // for MAX_LINE_LEN;
#include "openmc/error.h"
#include "openmc/tallies/filter_azimuthal.h"
#include "openmc/tallies/filter_cell.h"
#include "openmc/tallies/filter_cellborn.h"
#include "openmc/tallies/filter_cellfrom.h"
#include "openmc/tallies/filter_distribcell.h"
#include "openmc/tallies/filter_energyfunc.h"
#include "openmc/tallies/filter_legendre.h"
#include "openmc/tallies/filter_material.h"
#include "openmc/tallies/filter_mesh.h"
#include "openmc/tallies/filter_meshsurface.h"
#include "openmc/tallies/filter_mu.h"
#include "openmc/tallies/filter_polar.h"
#include "openmc/tallies/filter_sph_harm.h"
#include "openmc/tallies/filter_sptl_legendre.h"
#include "openmc/tallies/filter_surface.h"
#include "openmc/tallies/filter_universe.h"
#include "openmc/tallies/filter_zernike.h"


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

std::vector<FilterMatch> filter_matches;
std::vector<Filter*> tally_filters;

//==============================================================================
// Non-member functions
//==============================================================================

void
free_memory_tally_c()
{
  #pragma omp parallel
  {
    filter_matches.clear();
  }

  for (Filter* filt : tally_filters) {delete filt;}
  tally_filters.clear();
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  FilterMatch* filter_match_pointer(int indx)
  {return &filter_matches[indx];}

  void
  filter_match_bins_push_back(FilterMatch* match, int val)
  {match->bins_.push_back(val);}

  void
  filter_match_weights_push_back(FilterMatch* match, double val)
  {match->weights_.push_back(val);}

  void
  filter_match_bins_clear(FilterMatch* match)
  {match->bins_.clear();}

  void
  filter_match_weights_clear(FilterMatch* match)
  {match->weights_.clear();}

  int
  filter_match_bins_size(FilterMatch* match)
  {return match->bins_.size();}

  int
  filter_match_bins_data(FilterMatch* match, int indx)
  {return match->bins_[indx-1];}

  double
  filter_match_weights_data(FilterMatch* match, int indx)
  {return match->weights_[indx-1];}

  void
  filter_match_bins_set_data(FilterMatch* match, int indx, int val)
  {match->bins_[indx-1] = val;}

  Filter*
  allocate_filter(const char* type)
  {
    std::string type_ {type};
    if (type_ == "azimuthal") {
      tally_filters.push_back(new AzimuthalFilter());
    } else if (type_ == "cell") {
      tally_filters.push_back(new CellFilter());
    } else if (type_ == "cellborn") {
      tally_filters.push_back(new CellbornFilter());
    } else if (type_ == "cellfrom") {
      tally_filters.push_back(new CellFromFilter());
    } else if (type_ == "distribcell") {
      tally_filters.push_back(new DistribcellFilter());
    } else if (type_ == "energyfunction") {
      tally_filters.push_back(new EnergyFunctionFilter());
    } else if (type_ == "legendre") {
      tally_filters.push_back(new LegendreFilter());
    } else if (type_ == "material") {
      tally_filters.push_back(new MaterialFilter());
    } else if (type_ == "mesh") {
      tally_filters.push_back(new MeshFilter());
    } else if (type_ == "meshsurface") {
      tally_filters.push_back(new MeshSurfaceFilter());
    } else if (type_ == "mu") {
      tally_filters.push_back(new MuFilter());
    } else if (type_ == "polar") {
      tally_filters.push_back(new PolarFilter());
    } else if (type_ == "surface") {
      tally_filters.push_back(new SurfaceFilter());
    } else if (type_ == "spatiallegendre") {
      tally_filters.push_back(new SpatialLegendreFilter());
    } else if (type_ == "sphericalharmonics") {
      tally_filters.push_back(new SphericalHarmonicsFilter());
    } else if (type_ == "universe") {
      tally_filters.push_back(new UniverseFilter());
    } else if (type_ == "zernike") {
      tally_filters.push_back(new ZernikeFilter());
    } else if (type_ == "zernikeradial") {
      tally_filters.push_back(new ZernikeRadialFilter());
    } else {
      return nullptr;
    }
    return tally_filters.back();
  }

  void filter_from_xml(Filter* filt, pugi::xml_node* node)
  {filt->from_xml(*node);}

  void
  filter_get_all_bins(Filter* filt, Particle* p, int estimator,
                      FilterMatch* match)
  {
    filt->get_all_bins(p, estimator, *match);
  }

  void filter_to_statepoint(Filter* filt, hid_t group)
  {filt->to_statepoint(group);}

  void filter_text_label(Filter* filt, int bin, char* label)
  {
    std::string label_str = filt->text_label(bin);
    int i = 0;
    for (; i < label_str.size() && i < MAX_LINE_LEN; i++)
      label[i] = label_str[i];
    label[i] = '\0';
  }

  void filter_initialize(Filter* filt) {filt->initialize();}

  int filter_n_bins(Filter* filt) {return filt->n_bins_;}

  int mesh_filter_get_mesh(MeshFilter* filt) {return filt->mesh_;}

  int sphharm_filter_get_cosine(SphericalHarmonicsFilter* filt)
  {return static_cast<int>(filt->cosine_);}
}

} // namespace openmc
