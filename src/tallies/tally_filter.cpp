#include "openmc/tallies/tally_filter.h"

#include <string>

#include "openmc/constants.h"  // for MAX_LINE_LEN;
#include "openmc/tallies/tally_filter_azimuthal.h"
#include "openmc/tallies/tally_filter_cell.h"
#include "openmc/tallies/tally_filter_cellborn.h"
#include "openmc/tallies/tally_filter_cellfrom.h"
#include "openmc/tallies/tally_filter_distribcell.h"
#include "openmc/tallies/tally_filter_legendre.h"
#include "openmc/tallies/tally_filter_material.h"
#include "openmc/tallies/tally_filter_mesh.h"
#include "openmc/tallies/tally_filter_meshsurface.h"
#include "openmc/tallies/tally_filter_mu.h"
#include "openmc/tallies/tally_filter_polar.h"
#include "openmc/tallies/tally_filter_sph_harm.h"
#include "openmc/tallies/tally_filter_surface.h"
#include "openmc/tallies/tally_filter_universe.h"


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

std::vector<TallyFilterMatch> filter_matches;
std::vector<TallyFilter*> tally_filters;

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

  for (TallyFilter* filt : tally_filters) {delete filt;}
  tally_filters.clear();
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  TallyFilterMatch* filter_match_pointer(int indx)
  {return &filter_matches[indx];}

  void
  filter_match_bins_push_back(TallyFilterMatch* match, int val)
  {match->bins.push_back(val);}

  void
  filter_match_weights_push_back(TallyFilterMatch* match, double val)
  {match->weights.push_back(val);}

  void
  filter_match_bins_clear(TallyFilterMatch* match)
  {match->bins.clear();}

  void
  filter_match_weights_clear(TallyFilterMatch* match)
  {match->weights.clear();}

  int
  filter_match_bins_size(TallyFilterMatch* match)
  {return match->bins.size();}

  int
  filter_match_bins_data(TallyFilterMatch* match, int indx)
  {return match->bins.at(indx-1);}

  double
  filter_match_weights_data(TallyFilterMatch* match, int indx)
  {return match->weights.at(indx-1);}

  void
  filter_match_bins_set_data(TallyFilterMatch* match, int indx, int val)
  {match->bins.at(indx-1) = val;}

  TallyFilter*
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
    } else if (type_ == "sphericalharmonics") {
      tally_filters.push_back(new SphericalHarmonicsFilter());
    } else if (type_ == "universe") {
      tally_filters.push_back(new UniverseFilter());
    } else {
      return nullptr;
    }
    return tally_filters.back();
  }

  void filter_from_xml(TallyFilter* filt, pugi::xml_node* node)
  {filt->from_xml(*node);}

  void
  filter_get_all_bins(TallyFilter* filt, Particle* p, int estimator,
                      TallyFilterMatch* match)
  {
    filt->get_all_bins(p, estimator, *match);
  }

  void filter_to_statepoint(TallyFilter* filt, hid_t group)
  {filt->to_statepoint(group);}

  void filter_text_label(TallyFilter* filt, int bin, char* label)
  {
    std::string label_str = filt->text_label(bin);
    int i = 0;
    for (; i < label_str.size() && i < MAX_LINE_LEN; i++)
      label[i] = label_str[i];
    label[i] = '\0';
  }

  void filter_initialize(TallyFilter* filt) {filt->initialize();}

  int filter_n_bins(TallyFilter* filt) {return filt->n_bins_;}

  void
  cell_filter_get_bins(CellFilter* filt, int32_t** cells, int32_t* n)
  {
    *cells = filt->cells_.data();
    *n = filt->cells_.size();
  }

  int legendre_filter_get_order(LegendreFilter* filt)
  {return filt->order_;}

  void
  legendre_filter_set_order(LegendreFilter* filt, int order)
  {
    filt->order_ = order;
    filt->n_bins_ = order + 1;
  }

  void
  material_filter_get_bins(MaterialFilter* filt, int32_t** bins, int32_t* n)
  {
    *bins = filt->materials_.data();
    *n = filt->materials_.size();
  }

  void
  material_filter_set_bins(MaterialFilter* filt, int32_t n, int32_t* bins)
  {
    filt->materials_.clear();
    filt->materials_.resize(n);
    for (int i = 0; i < n; i++) filt->materials_[i] = bins[i];
    filt->n_bins_ = filt->materials_.size();
    filt->map_.clear();
    for (int i = 0; i < n; i++) filt->map_[filt->materials_[i]] = i;
  }

  int mesh_filter_get_mesh(MeshFilter* filt) {return filt->mesh_;}

  void
  mesh_filter_set_mesh(MeshFilter* filt, int mesh)
  {
    filt->mesh_ = mesh;
    filt->n_bins_ = 1;
    for (auto dim : meshes[mesh]->shape_) filt->n_bins_ *= dim;
  }

  int meshsurface_filter_get_mesh(MeshSurfaceFilter* filt) {return filt->mesh_;}

  void
  meshsurface_filter_set_mesh(MeshSurfaceFilter* filt, int mesh)
  {
    filt->mesh_ = mesh;
    filt->n_bins_ = 4 * meshes[mesh]->n_dimension_;
    for (auto dim : meshes[mesh]->shape_) filt->n_bins_ *= dim;
  }

  int sphharm_filter_get_order(SphericalHarmonicsFilter* filt)
  {return filt->order_;}

  int sphharm_filter_get_cosine(SphericalHarmonicsFilter* filt)
  {return static_cast<int>(filt->cosine_);}

  void
  sphharm_filter_set_order(SphericalHarmonicsFilter* filt, int order)
  {
    filt->order_ = order;
    filt->n_bins_ = (order + 1) * (order + 1);
  }

  void sphharm_filter_set_cosine(SphericalHarmonicsFilter* filt, int cosine)
  {filt->cosine_ = static_cast<SphericalHarmonicsCosine>(cosine);}
}

} // namespace openmc
