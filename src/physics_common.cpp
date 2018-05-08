#include "openmc/physics_common.h"
#include "openmc/constants.h"
#include "openmc/settings.h"
#include "openmc/random_lcg.h"
#include "openmc/variance_reduction.h"

namespace openmc {

//==============================================================================
// RUSSIAN_ROULETTE
//==============================================================================

void russian_roulette(Particle* p)
{
  if (p->wgt_ < settings::weight_cutoff) {
    if (prn() < p->wgt_ / settings::weight_survive) {
      p->wgt_ = settings::weight_survive;
      p->wgt_last_ = p->wgt_;
    } else {
      p->wgt_ = 0.;
      p->wgt_last_ = 0.;
      p->alive_ = false;
    }
  }
}

//==============================================================================
// IMPORTANCE_SPLIT
//==============================================================================

void split_by_ratio(Particle* p, double importance_ratio) {
  p->last_wgt = p->wgt;
  double importance_rounded = std::floor(importance_ratio);
  
  int num_splits = 0;
  // how many splits - probabilistic 3/2 split 
  //  half the time we get 3 and half we get 2
  if ( prn() < importance_ratio - importance_rounded )
    num_splits = int(importance_rounded);
  else 
    num_splits = int(std::ceil(importance_ratio));
  // calculate the tru number of splits
  p->wgt /= float(num_splits);
  // note 1 -> n not 0->n
  for ( int i = 1 ; i < num_splits ; i++ ) {
    p->create_secondary(p->uvw,p->E,p->type,true);
  }

  return;
}

//==============================================================================
// ROULETTE
//==============================================================================

void roulette(Particle* p, double importance_ratio) {
  // this is a kill event
  if ( prn() >= importance_ratio) {
    p->alive = false;
    p->wgt = 0.;
    p->last_wgt = 0.;
  } else {
    p->wgt /= importance_ratio
  }
  return;
}

//==============================================================================
// WEIGHT_WINDOW - assume we have already checked cutoffs
//==============================================================================

void weight_window_split(Particle *p)

  if ( mesh_weight) 
    get_window_from_mesh(p,weight_low,weight_weight_high);
  else
    get_window_from_cell(p,weight_low,weight_weight_high)
  
  if ( p->wgt > weight_high) 
    if ( p->wgt/weight_high > MAX_SPLIT ) {
          // long history mitigation - if the number of requested
          // splits is greater than the prescribed amount in the input
          // lets say a million then we increase the effective weight
          // map of the problem until this history ends

          // todo - dynamic sharing of the split history oversplitting
          // is an issue mostly because of the parallel efficiency drop
          // perhaps this seconday particle can be distributed over the
          // source bank instead? need to be sure the random number state
          // isnt peturbed      
    }
    /* 
    // split into n secondaries
    int n_secondary = std::floor(p->wgt/weight_high);
    double secondary_wgt = p->wgt/n_secondary;
    double remaining_wgt = p->wgt - n_secondary*secondary_wgt;
  
    p->wgt = secondary_wgt;
    // create secondaries with the correct weight
    for ( int i = 1 ; i < n_secondary) {
      p->create_secondary(p->uvw,p->E,p->type,true);
    }
    // reset original particle to include left over weight
    p->wgt += remaining_wgt;
    */
    split_by_ratio(p,p->wgt/weight_high);
  else if ( p->wgt < weight_low ) {
    roulette(p,p->wgt/weight_low)  
  } else {
    // do nothing we are in the window
  }
  return;  
}

// split according to importance 
void importance_split(Particle *p){

  // importances can only come from cells
  int current_cell = p->last_n_coord[p->n_coord]->cell;
  int last_cell = p->last_n_coord[p->n_coord]->last_cell;

  double importance_ratio = variance_reduction::importances[current_cell]/
                            variance_reduciton::importances[last_cell];

  if (importance_ratio > 1) 
    split_by_ratio(p, importance_ratio);
  else
    roulette(p, importance_ratio);  
}


//==============================================================================
// VARIANCE_REDUCTION
//==============================================================================

void variance_reduction(Particle* p)
{
  // check weight cutoff
  if ( p->wgt < settings::weight_cutoff ) {
    p->wgt = 0.;
    p->last_wgt = 0.;
    p->alive = false;
    return;
  }

  if ( importance_splitting ) 
    importance_split(p);
  else 
    weight_window_split(p);
}

} //namespace openmc
