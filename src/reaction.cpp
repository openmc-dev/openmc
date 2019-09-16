#include "openmc/reaction.h"

#include <string>
#include <utility> // for move

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/endf.h"
#include "openmc/random_lcg.h"
#include "openmc/secondary_uncorrelated.h"

namespace openmc {

//==============================================================================
// Reaction implementation
//==============================================================================

Reaction::Reaction(hid_t group, const std::vector<int>& temperatures)
{
  read_attribute(group, "Q_value", q_value_);
  read_attribute(group, "mt", mt_);
  int tmp;
  read_attribute(group, "center_of_mass", tmp);
  scatter_in_cm_ = (tmp == 1);

  // Checks if redudant attribute exists before loading
  // (for compatibiltiy with legacy .h5 libraries)
  if (attribute_exists(group, "redundant")) {
    read_attribute(group, "redundant", tmp);
    redundant_ = (tmp == 1);
  } else {
    redundant_ = false;
  }

  // Read cross section and threshold_idx data
  for (auto t : temperatures) {
    // Get group corresponding to temperature
    std::string temp_str {std::to_string(t) + "K"};
    hid_t temp_group = open_group(group, temp_str.c_str());
    hid_t dset = open_dataset(temp_group, "xs");

    // Get threshold index
    TemperatureXS xs;
    read_attribute(dset, "threshold_idx", xs.threshold);

    // Read cross section values
    read_dataset(dset, xs.value);
    close_dataset(dset);
    close_group(temp_group);

    // create new entry in xs vector
    xs_.push_back(std::move(xs));
  }

  // Read products
  for (const auto& name : group_names(group)) {
    if (name.rfind("product_", 0) == 0) {
      hid_t pgroup = open_group(group, name.c_str());
      products_.emplace_back(pgroup);
      close_group(pgroup);
    }
  }

  // <<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<
  // Before the secondary distribution refactor, when the angle/energy
  // distribution was uncorrelated, no angle was actually sampled. With
  // the refactor, an angle is always sampled for an uncorrelated
  // distribution even when no angle distribution exists in the ACE file
  // (isotropic is assumed). To preserve the RNG stream, we explicitly
  // mark fission reactions so that we avoid the angle sampling.
  if (is_fission(mt_)) {
    for (auto& p : products_) {
      if (p.particle_ == Particle::Type::neutron) {
        for (auto& d : p.distribution_) {
          auto d_ = dynamic_cast<UncorrelatedAngleEnergy*>(d.get());
          if (d_) d_->fission() = true;
        }
      }
    }
  }
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<
}

//==============================================================================
// Non-member functions
//==============================================================================

std::string reaction_name(int mt)
{
  if (mt == SCORE_FLUX) {
    return "flux";
  } else if (mt == SCORE_TOTAL) {
    return "total";
  } else if (mt == SCORE_SCATTER) {
    return "scatter";
  } else if (mt == SCORE_NU_SCATTER) {
    return "nu-scatter";
  } else if (mt == SCORE_ABSORPTION) {
    return "absorption";
  } else if (mt == SCORE_FISSION) {
    return "fission";
  } else if (mt == SCORE_NU_FISSION) {
    return "nu-fission";
  } else if (mt == SCORE_DECAY_RATE) {
    return "decay-rate";
  } else if (mt == SCORE_DELAYED_NU_FISSION) {
    return "delayed-nu-fission";
  } else if (mt == SCORE_PROMPT_NU_FISSION) {
    return "prompt-nu-fission";
  } else if (mt == SCORE_KAPPA_FISSION) {
    return "kappa-fission";
  } else if (mt == SCORE_CURRENT) {
    return "current";
  } else if (mt == SCORE_EVENTS) {
    return "events";
  } else if (mt == SCORE_INVERSE_VELOCITY) {
    return "inverse-velocity";
  } else if (mt == SCORE_FISS_Q_PROMPT) {
    return "fission-q-prompt";
  } else if (mt == SCORE_FISS_Q_RECOV) {
    return "fission-q-recoverable";

  // Normal ENDF-based reactions
  } else if (mt == TOTAL_XS) {
    return "(n,total)";
  } else if (mt == ELASTIC) {
    return "(n,elastic)";
  } else if (mt == N_LEVEL) {
    return "(n,level)";
  } else if (mt == N_2ND) {
    return "(n,2nd)";
  } else if (mt == N_2N) {
    return "(n,2n)";
  } else if (mt == N_3N) {
    return "(n,3n)";
  } else if (mt == N_FISSION) {
    return "(n,fission)";
  } else if (mt == N_F) {
    return "(n,f)";
  } else if (mt == N_NF) {
    return "(n,nf)";
  } else if (mt == N_2NF) {
    return "(n,2nf)";
  } else if (mt == N_NA) {
    return "(n,na)";
  } else if (mt == N_N3A) {
    return "(n,n3a)";
  } else if (mt == N_2NA) {
    return "(n,2na)";
  } else if (mt == N_3NA) {
    return "(n,3na)";
  } else if (mt == N_NP) {
    return "(n,np)";
  } else if (mt == N_N2A) {
    return "(n,n2a)";
  } else if (mt == N_2N2A) {
    return "(n,2n2a)";
  } else if (mt == N_ND) {
    return "(n,nd)";
  } else if (mt == N_NT) {
    return "(n,nt)";
  } else if (mt == N_N3HE) {
    return "(n,nHe-3)";
  } else if (mt == N_ND2A) {
    return "(n,nd2a)";
  } else if (mt == N_NT2A) {
    return "(n,nt2a)";
  } else if (mt == N_4N) {
    return "(n,4n)";
  } else if (mt == N_3NF) {
    return "(n,3nf)";
  } else if (mt == N_2NP) {
    return "(n,2np)";
  } else if (mt == N_3NP) {
    return "(n,3np)";
  } else if (mt == N_N2P) {
    return "(n,n2p)";
  } else if (mt == N_NPA) {
    return "(n,npa)";
  } else if (N_N1 <= mt && mt <= N_N40) {
    return "(n,n" + std::to_string(mt-50) + ")";
  } else if (mt == N_NC) {
    return "(n,nc)";
  } else if (mt == N_DISAPPEAR) {
    return "(n,disappear)";
  } else if (mt == N_GAMMA) {
    return "(n,gamma)";
  } else if (mt == N_P) {
    return "(n,p)";
  } else if (mt == N_D) {
    return "(n,d)";
  } else if (mt == N_T) {
    return "(n,t)";
  } else if (mt == N_3HE) {
    return "(n,3He)";
  } else if (mt == N_A) {
    return "(n,a)";
  } else if (mt == N_2A) {
    return "(n,2a)";
  } else if (mt == N_3A) {
    return "(n,3a)";
  } else if (mt == N_2P) {
    return "(n,2p)";
  } else if (mt == N_PA) {
    return "(n,pa)";
  } else if (mt == N_T2A) {
    return "(n,t2a)";
  } else if (mt == N_D2A) {
    return "(n,d2a)";
  } else if (mt == N_PD) {
    return "(n,pd)";
  } else if (mt == N_PT) {
    return "(n,pt)";
  } else if (mt == N_DA) {
    return "(n,da)";
  } else if (mt == 201) {
    return "(n,Xn)";
  } else if (mt == 202) {
    return "(n,Xgamma)";
  } else if (mt == N_XP) {
    return "(n,Xp)";
  } else if (mt == N_XD) {
    return "(n,Xd)";
  } else if (mt == N_XT) {
    return "(n,Xt)";
  } else if (mt == N_X3HE) {
    return "(n,X3He)";
  } else if (mt == N_XA) {
    return "(n,Xa)";
  } else if (mt == HEATING) {
    return "heating";
  } else if (mt == DAMAGE_ENERGY) {
    return "damage-energy";
  } else if (mt == COHERENT) {
    return "coherent scatter";
  } else if (mt == INCOHERENT) {
    return "incoherent scatter";
  } else if (mt == PAIR_PROD_ELEC) {
    return "pair production, electron";
  } else if (mt == PAIR_PROD) {
    return "pair production";
  } else if (mt == PAIR_PROD_NUC) {
    return "pair production, nuclear";
  } else if (mt == PHOTOELECTRIC) {
    return "photoelectric";
  } else if (534 <= mt && mt <= 572) {
    std::stringstream name;
    name << "photoelectric, " << SUBSHELLS[mt - 534] << " subshell";
    return name.str();
  } else if (600 <= mt && mt <= 648) {
    return "(n,p" + std::to_string(mt-600) + ")";
  } else if (mt == 649) {
    return "(n,pc)";
  } else if (650 <= mt && mt <= 698) {
    return "(n,d" + std::to_string(mt-650) + ")";
  } else if (mt == 699) {
    return "(n,dc)";
  } else if (700 <= mt && mt <= 748) {
    return "(n,t" + std::to_string(mt-700) + ")";
  } else if (mt == 749) {
    return "(n,tc)";
  } else if (750 <= mt && mt <= 798) {
    return "(n,3He" + std::to_string(mt-750) + ")";
  } else if (mt == 799) {
    return "(n,3Hec)";
  } else if (800 <= mt && mt <= 848) {
    return "(n,a" + std::to_string(mt-800) + ")";
  } else if (mt == 849) {
    return "(n,ac)";
  } else if (mt == HEATING_LOCAL) {
    return "heating-local";
  } else {
    return "MT=" + std::to_string(mt);
  }
}

} // namespace openmc
