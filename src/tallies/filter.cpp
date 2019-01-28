#include "openmc/tallies/filter.h"

#include <string>

#include "openmc/capi.h"
#include "openmc/constants.h"  // for MAX_LINE_LEN;
#include "openmc/error.h"
#include "openmc/tallies/filter_azimuthal.h"
#include "openmc/tallies/filter_cell.h"
#include "openmc/tallies/filter_cellborn.h"
#include "openmc/tallies/filter_cellfrom.h"
#include "openmc/tallies/filter_delayedgroup.h"
#include "openmc/tallies/filter_distribcell.h"
#include "openmc/tallies/filter_energyfunc.h"
#include "openmc/tallies/filter_energy.h"
#include "openmc/tallies/filter_legendre.h"
#include "openmc/tallies/filter_material.h"
#include "openmc/tallies/filter_mesh.h"
#include "openmc/tallies/filter_meshsurface.h"
#include "openmc/tallies/filter_mu.h"
#include "openmc/tallies/filter_particle.h"
#include "openmc/tallies/filter_polar.h"
#include "openmc/tallies/filter_sph_harm.h"
#include "openmc/tallies/filter_sptl_legendre.h"
#include "openmc/tallies/filter_surface.h"
#include "openmc/tallies/filter_universe.h"
#include "openmc/tallies/filter_zernike.h"

// explicit template instantiation definition
template class std::vector<openmc::FilterMatch>;

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {
std::vector<FilterMatch> filter_matches;
} // namespace simulation

namespace model {
std::vector<std::unique_ptr<Filter>> tally_filters;
} // namespace model

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  // filter_match_point moved to simulation.cpp

  int
  filter_match_get_i_bin(FilterMatch* match)
  {return match->i_bin_;}

  void
  filter_match_set_i_bin(FilterMatch* match, int i)
  {match->i_bin_ = i;}

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
      model::tally_filters.push_back(std::make_unique<AzimuthalFilter>());
    } else if (type_ == "cell") {
      model::tally_filters.push_back(std::make_unique<CellFilter>());
    } else if (type_ == "cellborn") {
      model::tally_filters.push_back(std::make_unique<CellbornFilter>());
    } else if (type_ == "cellfrom") {
      model::tally_filters.push_back(std::make_unique<CellFromFilter>());
    } else if (type_ == "distribcell") {
      model::tally_filters.push_back(std::make_unique<DistribcellFilter>());
    } else if (type_ == "delayedgroup") {
      model::tally_filters.push_back(std::make_unique<DelayedGroupFilter>());
    } else if (type_ == "energyfunction") {
      model::tally_filters.push_back(std::make_unique<EnergyFunctionFilter>());
    } else if (type_ == "energy") {
      model::tally_filters.push_back(std::make_unique<EnergyFilter>());
    } else if (type_ == "energyout") {
      model::tally_filters.push_back(std::make_unique<EnergyoutFilter>());
    } else if (type_ == "legendre") {
      model::tally_filters.push_back(std::make_unique<LegendreFilter>());
    } else if (type_ == "material") {
      model::tally_filters.push_back(std::make_unique<MaterialFilter>());
    } else if (type_ == "mesh") {
      model::tally_filters.push_back(std::make_unique<MeshFilter>());
    } else if (type_ == "meshsurface") {
      model::tally_filters.push_back(std::make_unique<MeshSurfaceFilter>());
    } else if (type_ == "mu") {
      model::tally_filters.push_back(std::make_unique<MuFilter>());
    } else if (type_ == "particle") {
      model::tally_filters.push_back(std::make_unique<ParticleFilter>());
    } else if (type_ == "polar") {
      model::tally_filters.push_back(std::make_unique<PolarFilter>());
    } else if (type_ == "surface") {
      model::tally_filters.push_back(std::make_unique<SurfaceFilter>());
    } else if (type_ == "spatiallegendre") {
      model::tally_filters.push_back(std::make_unique<SpatialLegendreFilter>());
    } else if (type_ == "sphericalharmonics") {
      model::tally_filters.push_back(std::make_unique<SphericalHarmonicsFilter>());
    } else if (type_ == "universe") {
      model::tally_filters.push_back(std::make_unique<UniverseFilter>());
    } else if (type_ == "zernike") {
      model::tally_filters.push_back(std::make_unique<ZernikeFilter>());
    } else if (type_ == "zernikeradial") {
      model::tally_filters.push_back(std::make_unique<ZernikeRadialFilter>());
    } else {
      return nullptr;
    }
    return model::tally_filters.back().get();
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

  int mesh_filter_get_mesh(MeshFilter* filt) {return filt->mesh();}

  int sphharm_filter_get_cosine(SphericalHarmonicsFilter* filt)
  {return static_cast<int>(filt->cosine_);}
}

} // namespace openmc
