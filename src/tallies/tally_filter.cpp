#include "openmc/tallies/tally_filter.h"

#include <string>

#include "openmc/tallies/tally_filter_cell.h"
#include "openmc/tallies/tally_filter_cellfrom.h"
#include "openmc/tallies/tally_filter_mesh.h"
#include "openmc/tallies/tally_filter_meshsurface.h"


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
    if (type_ == "cell") {
      tally_filters.push_back(new CellFilter());
    } else if (type_ == "cellfrom") {
      tally_filters.push_back(new CellFromFilter());
    } else if (type_ == "mesh") {
      tally_filters.push_back(new MeshFilter());
    } else if (type_ == "meshsurface") {
      tally_filters.push_back(new MeshSurfaceFilter());
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

  void filter_initialize(TallyFilter* filt) {filt->initialize();}
}

} // namespace openmc
