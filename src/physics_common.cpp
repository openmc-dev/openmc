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
  p->wgt_last_ = p->wgt_;
  double importance_rounded = std::floor(importance_ratio);
  
  int num_splits = 0;
  // how many splits - probabilistic 3/2 split 
  //  half the time we get 3 and half we get 2
  if ( prn() < importance_ratio - importance_rounded )
    num_splits = int(importance_rounded);
  else 
    num_splits = int(std::ceil(importance_ratio));
  // calculate the tru number of splits
  p->wgt_ /= float(num_splits);
  // note 1 -> n not 0->n
  for ( int i = 1 ; i < num_splits ; i++ ) {
    p->create_secondary(p->u(),p->E_,p->type_);
  }

  return;
}

//==============================================================================
// ROULETTE
//==============================================================================

void roulette(Particle* p, double importance_ratio) {
  // this is a kill event
  double random = prn();
  //std::cout << random << " " << importance_ratio;
  if ( random < importance_ratio) {
    //std::cout << " survive " << std::endl;
    p->wgt_last_ = p->wgt_;
    //std::cout << p->last_wgt << " ";
    p->wgt_ /= importance_ratio;
    //std::cout << p->wgt << std::endl;  
  } else {
    //std::cout << " kill" << std::endl;
    p->alive_ = false;
    p->wgt_ = 0.;
    p->wgt_last_ = 0.;
  }
  return;
}

//==============================================================================
// WEIGHT_WINDOW - assume we have already checked cutoffs
//==============================================================================

void weight_window_split(Particle *p) {
  bool mesh_weight = false;
  double weight_low = 0.;
  double weight_high = 0.;

  /* - todo complete weight from mesh
  if ( mesh_weight ) 
    get_window_from_mesh(p,weight_low,weight_high);
  else
    get_window_from_cell(p,weight_low,weight_high);
  */

  if ( p->wgt_ > weight_high ) {
    if ( p->wgt_/weight_high > MAX_SPLIT ) {
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
    split_by_ratio(p,p->wgt_/weight_high);
  } else if ( p->wgt_ < weight_low ) {
    roulette(p,p->wgt_/weight_low);  
  } else {
    // do nothing we are in the window
  }
  return;  
}

// split according to importance 
void importance_split(Particle *p){

  // importances can only come from cells
  int current_cell = p->coord_[p->n_coord_-1].cell;
  int last_cell = p->cell_last_[p->n_coord_last_-1];

  // if we havent changed cell
  //std::cout << p->id << " " << current_cell << " " << last_cell << std::endl;

  if ( ((current_cell < 0 ) ? (0) : (current_cell)) == 
       ((last_cell < 0) ? (0) : (last_cell))) return;

  //std::cout << p->id << " " << current_cell << " " << last_cell << std::endl;

  double importance_ratio = variance_reduction::importances[current_cell+1]/
                            variance_reduction::importances[last_cell+1];
  //std::cout << variance_reduction::importances[current_cell+1] << " ";
  //std::cout << variance_reduction::importances[last_cell+1] << std::endl;
  //std::cout << importance_ratio << std::endl;
  if (importance_ratio > 1) {
    split_by_ratio(p, importance_ratio);
  } else if ( importance_ratio < 1 ) {
    roulette(p, importance_ratio);
  }  
  return;
}


//==============================================================================
// VARIANCE_REDUCTION
//==============================================================================

void perform_vr(Particle* p)
{
  // check weight cutoff
  if ( p->wgt_ < settings::weight_cutoff ) {
    p->wgt_ = 0.;
    p->wgt_last_ = 0.;
    p->alive_ = false;
    return;
  }

  if (p->wgt_absorb_ == 0) return;

  if (variance_reduction::importance_splitting) 
    importance_split(p);
  else if (variance_reduction::weight_splitting)
    weight_window_split(p);
}

} //namespace openmc
