#include "openmc/tallies/tally_filter.h"

#include <string>

#include "openmc/capi.h"
#include "openmc/constants.h"  // for MAX_LINE_LEN;
#include "openmc/error.h"
#include "openmc/tallies/tally_filter_azimuthal.h"
#include "openmc/tallies/tally_filter_cell.h"
#include "openmc/tallies/tally_filter_cellborn.h"
#include "openmc/tallies/tally_filter_cellfrom.h"
#include "openmc/tallies/tally_filter_distribcell.h"
#include "openmc/tallies/tally_filter_energyfunc.h"
#include "openmc/tallies/tally_filter_legendre.h"
#include "openmc/tallies/tally_filter_material.h"
#include "openmc/tallies/tally_filter_mesh.h"
#include "openmc/tallies/tally_filter_meshsurface.h"
#include "openmc/tallies/tally_filter_mu.h"
#include "openmc/tallies/tally_filter_polar.h"
#include "openmc/tallies/tally_filter_sph_harm.h"
#include "openmc/tallies/tally_filter_sptl_legendre.h"
#include "openmc/tallies/tally_filter_surface.h"
#include "openmc/tallies/tally_filter_universe.h"
#include "openmc/tallies/tally_filter_zernike.h"


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
// C-API functions
//==============================================================================

// Fortran functions that will be called from C++
extern "C" int verify_filter(int32_t index);
extern "C" TallyFilter* filter_from_f(int32_t index);
extern "C" void filter_update_n_bins(int32_t index);

extern "C" {
  int
  openmc_cell_filter_get_bins(int32_t index, int32_t** cells, int32_t* n)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "cell") {
      set_errmsg("Tried to get cells from a non-cell filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto cell_filt = static_cast<CellFilter*>(filt);
    *cells = cell_filt->cells_.data();
    *n = cell_filt->cells_.size();
    return 0;
  }

  int
  openmc_legendre_filter_get_order(int32_t index, int* order)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "legendre") {
      set_errmsg("Tried to get order on a non-expansion filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto l_filt = static_cast<LegendreFilter*>(filt);
    *order = l_filt->order_;
    return 0;
  }

  int
  openmc_legendre_filter_set_order(int32_t index, int order)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "legendre") {
      set_errmsg("Tried to set order on a non-expansion filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto l_filt = static_cast<LegendreFilter*>(filt);
    l_filt->order_ = order;
    l_filt->n_bins_ = order + 1;
    filter_update_n_bins(index);
    return 0;
  }

  int
  openmc_material_filter_get_bins(int32_t index, int32_t** bins, int32_t* n)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "material") {
      set_errmsg("Tried to get material filter bins on a non-material filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto mat_filt = static_cast<MaterialFilter*>(filt);
    *bins = mat_filt->materials_.data();
    *n = mat_filt->materials_.size();
    return 0;
  }

  int
  openmc_material_filter_set_bins(int32_t index, int32_t n, const int32_t* bins)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "material") {
      set_errmsg("Tried to set material filter bins on a non-material filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto mat_filt = static_cast<MaterialFilter*>(filt);
    mat_filt->materials_.clear();
    mat_filt->materials_.resize(n);
    for (int i = 0; i < n; i++) mat_filt->materials_[i] = bins[i];
    mat_filt->n_bins_ = mat_filt->materials_.size();
    mat_filt->map_.clear();
    for (int i = 0; i < n; i++) mat_filt->map_[mat_filt->materials_[i]] = i;
    filter_update_n_bins(index);
    return 0;
  }

  int
  openmc_mesh_filter_get_mesh(int32_t index, int32_t* index_mesh)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "mesh" && filt->type() != "meshsurface") {
      set_errmsg("Tried to get mesh on a non-mesh filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto mesh_filt = static_cast<MeshFilter*>(filt);
    *index_mesh = mesh_filt->mesh_;
    return 0;
  }

  int
  openmc_mesh_filter_set_mesh(int32_t index, int32_t index_mesh)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "mesh" && filt->type() != "meshsurface") {
      set_errmsg("Tried to set mesh on a non-mesh filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    if (index_mesh < 0 || index_mesh >= meshes.size()) {
      set_errmsg("Index in 'meshes' array is out of bounds.");
      return OPENMC_E_OUT_OF_BOUNDS;
    }

    auto mesh_filt = static_cast<MeshFilter*>(filt);
    mesh_filt->mesh_ = index_mesh;
    if (filt->type() == "mesh") {
      mesh_filt->n_bins_ = 1;
    } else {
      filt->n_bins_ = 4 * meshes[index_mesh]->n_dimension_;
    }
    for (auto dim : meshes[index_mesh]->shape_) mesh_filt->n_bins_ *= dim;
    filter_update_n_bins(index);
    return 0;
  }

  int openmc_meshsurface_filter_get_mesh(int32_t index, int32_t* index_mesh)
  {return openmc_mesh_filter_get_mesh(index, index_mesh);}

  int openmc_meshsurface_filter_set_mesh(int32_t index, int32_t index_mesh)
  {return openmc_mesh_filter_set_mesh(index, index_mesh);}

  int
  openmc_spatial_legendre_filter_get_order(int32_t index, int* order)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "spatiallegendre") {
      set_errmsg("Not a spatial Legendre filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto l_filt = static_cast<SpatialLegendreFilter*>(filt);
    *order = l_filt->order_;
    return 0;
  }

  int
  openmc_spatial_legendre_filter_get_params(int32_t index, int* axis,
                                            double* min, double* max)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "spatiallegendre") {
      set_errmsg("Not a spatial Legendre filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto l_filt = static_cast<SpatialLegendreFilter*>(filt);
    *axis = static_cast<int>(l_filt->axis_);
    *min = l_filt->min_;
    *max = l_filt->max_;
    return 0;
  }

  int
  openmc_spatial_legendre_filter_set_order(int32_t index, int order)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "spatiallegendre") {
      set_errmsg("Not a spatial Legendre filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto l_filt = static_cast<SpatialLegendreFilter*>(filt);
    l_filt->order_ = order;
    l_filt->n_bins_ = order + 1;
    filter_update_n_bins(index);
    return 0;
  }

  int
  openmc_spatial_legendre_filter_set_params(int32_t index, const int* axis,
    const double* min, const double* max)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "spatiallegendre") {
      set_errmsg("Not a spatial Legendre filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto l_filt = static_cast<SpatialLegendreFilter*>(filt);
    if (axis) l_filt->axis_ = static_cast<LegendreAxis>(*axis);
    if (min) l_filt->min_ = *min;
    if (max) l_filt->max_ = *max;
    return 0;
  }

  int
  openmc_sphharm_filter_get_order(int32_t index, int* order)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "sphericalharmonics") {
      set_errmsg("Not a spherical harmonics filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto sph_filt = static_cast<SphericalHarmonicsFilter*>(filt);
    *order = sph_filt->order_;
    return 0;
  }

  int
  openmc_sphharm_filter_get_cosine(int32_t index, char cosine[])
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "sphericalharmonics") {
      set_errmsg("Not a spherical harmonics filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto sph_filt = static_cast<SphericalHarmonicsFilter*>(filt);
    if (sph_filt->cosine_ == SphericalHarmonicsCosine::scatter) {
      strcpy(cosine, "scatter");
    } else {
      strcpy(cosine, "particle");
    }
    return 0;
  }

  int
  openmc_sphharm_filter_set_order(int32_t index, int order)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "sphericalharmonics") {
      set_errmsg("Not a spherical harmonics filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto sph_filt = static_cast<SphericalHarmonicsFilter*>(filt);
    sph_filt->order_ = order;
    sph_filt->n_bins_ = (order + 1) * (order + 1);
    filter_update_n_bins(index);
    return 0;
  }

  int
  openmc_sphharm_filter_set_cosine(int32_t index, const char cosine[])
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "sphericalharmonics") {
      set_errmsg("Not a spherical harmonics filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto sph_filt = static_cast<SphericalHarmonicsFilter*>(filt);
    if (strcmp(cosine, "scatter") == 0) {
      sph_filt->cosine_ = SphericalHarmonicsCosine::scatter;
    } else if (strcmp(cosine, "particle") == 0) {
      sph_filt->cosine_ = SphericalHarmonicsCosine::particle;
    } else {
      set_errmsg("Invalid spherical harmonics cosine.");
      return OPENMC_E_INVALID_ARGUMENT;
    }
    return 0;
  }

  int
  openmc_zernike_filter_get_order(int32_t index, int* order)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "zernike" && filt->type() != "zernikeradial") {
      set_errmsg("Not a Zernike filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto z_filt = static_cast<ZernikeFilter*>(filt);
    *order = z_filt->order_;
    return 0;
  }

  int
  openmc_zernike_filter_get_params(int32_t index, double* x, double* y,
                                   double* r)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "zernike" && filt->type() != "zernikeradial") {
      set_errmsg("Not a Zernike filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto z_filt = static_cast<ZernikeFilter*>(filt);
    *x = z_filt->x_;
    *y = z_filt->y_;
    *r = z_filt->r_;
    return 0;
  }

  int
  openmc_zernike_filter_set_order(int32_t index, int order)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "zernike" && filt->type() != "zernikeradial") {
      set_errmsg("Not a Zernike filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto z_filt = static_cast<ZernikeFilter*>(filt);
    z_filt->order_ = order;
    z_filt->calc_n_bins();
    filter_update_n_bins(index);
    return 0;
  }

  int
  openmc_zernike_filter_set_params(int32_t index, const double* x,
                                   const double* y, const double* r)
  {
    int err = verify_filter(index);
    if (err) return err;

    auto filt = filter_from_f(index);
    if (filt->type() != "zernike" && filt->type() != "zernikeradial") {
      set_errmsg("Not a Zernike filter.");
      return OPENMC_E_INVALID_TYPE;
    }

    auto z_filt = static_cast<ZernikeFilter*>(filt);
    if (x) z_filt->x_ = *x;
    if (y) z_filt->y_ = *y;
    if (r) z_filt->r_ = *r;
    return 0;
  }
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  TallyFilterMatch* filter_match_pointer(int indx)
  {return &filter_matches[indx];}

  void
  filter_match_bins_push_back(TallyFilterMatch* match, int val)
  {match->bins_.push_back(val);}

  void
  filter_match_weights_push_back(TallyFilterMatch* match, double val)
  {match->weights_.push_back(val);}

  void
  filter_match_bins_clear(TallyFilterMatch* match)
  {match->bins_.clear();}

  void
  filter_match_weights_clear(TallyFilterMatch* match)
  {match->weights_.clear();}

  int
  filter_match_bins_size(TallyFilterMatch* match)
  {return match->bins_.size();}

  int
  filter_match_bins_data(TallyFilterMatch* match, int indx)
  {return match->bins_[indx-1];}

  double
  filter_match_weights_data(TallyFilterMatch* match, int indx)
  {return match->weights_[indx-1];}

  void
  filter_match_bins_set_data(TallyFilterMatch* match, int indx, int val)
  {match->bins_[indx-1] = val;}

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

  int mesh_filter_get_mesh(MeshFilter* filt) {return filt->mesh_;}

  int sphharm_filter_get_cosine(SphericalHarmonicsFilter* filt)
  {return static_cast<int>(filt->cosine_);}
}

} // namespace openmc
