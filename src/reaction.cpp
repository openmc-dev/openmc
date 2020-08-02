#include "openmc/reaction.h"

#include <string>
#include <unordered_map>
#include <utility> // for move

#include <fmt/core.h>

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
    hid_t temp_group = open_group(group, fmt::format("{}K", t).c_str());
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

const std::unordered_map<int, std::string> REACTION_NAME_MAP {
  {SCORE_FLUX, "flux"},
  {SCORE_TOTAL, "total"},
  {SCORE_SCATTER, "scatter"},
  {SCORE_NU_SCATTER, "nu-scatter"},
  {SCORE_ABSORPTION, "absorption"},
  {SCORE_FISSION, "fission"},
  {SCORE_NU_FISSION, "nu-fission"},
  {SCORE_DECAY_RATE, "decay-rate"},
  {SCORE_DELAYED_NU_FISSION, "delayed-nu-fission"},
  {SCORE_PROMPT_NU_FISSION, "prompt-nu-fission"},
  {SCORE_KAPPA_FISSION, "kappa-fission"},
  {SCORE_CURRENT, "current"},
  {SCORE_EVENTS, "events"},
  {SCORE_INVERSE_VELOCITY, "inverse-velocity"},
  {SCORE_FISS_Q_PROMPT, "fission-q-prompt"},
  {SCORE_FISS_Q_RECOV, "fission-q-recoverable"},
  // Normal ENDF-based reactions
  {TOTAL_XS, "(n,total)"},
  {ELASTIC, "(n,elastic)"},
  {N_LEVEL, "(n,level)"},
  {N_2ND, "(n,2nd)"},
  {N_2N, "(n,2n)"},
  {N_3N, "(n,3n)"},
  {N_FISSION, "(n,fission)"},
  {N_F, "(n,f)"},
  {N_NF, "(n,nf)"},
  {N_2NF, "(n,2nf)"},
  {N_NA, "(n,na)"},
  {N_N3A, "(n,n3a)"},
  {N_2NA, "(n,2na)"},
  {N_3NA, "(n,3na)"},
  {N_NP, "(n,np)"},
  {N_N2A, "(n,n2a)"},
  {N_2N2A, "(n,2n2a)"},
  {N_ND, "(n,nd)"},
  {N_NT, "(n,nt)"},
  {N_N3HE, "(n,n3He)"},
  {N_ND2A, "(n,nd2a)"},
  {N_NT2A, "(n,nt2a)"},
  {N_4N, "(n,4n)"},
  {N_3NF, "(n,3nf)"},
  {N_2NP, "(n,2np)"},
  {N_3NP, "(n,3np)"},
  {N_N2P, "(n,n2p)"},
  {N_NPA, "(n,npa)"},
  {N_NC, "(n,nc)"},
  {N_DISAPPEAR, "(n,disappear)"},
  {N_GAMMA, "(n,gamma)"},
  {N_P, "(n,p)"},
  {N_D, "(n,d)"},
  {N_T, "(n,t)"},
  {N_3HE, "(n,3He)"},
  {N_A, "(n,a)"},
  {N_2A, "(n,2a)"},
  {N_3A, "(n,3a)"},
  {N_2P, "(n,2p)"},
  {N_PA, "(n,pa)"},
  {N_T2A, "(n,t2a)"},
  {N_D2A, "(n,d2a)"},
  {N_PD, "(n,pd)"},
  {N_PT, "(n,pt)"},
  {N_DA, "(n,da)"},
  {201, "(n,Xn)"},
  {202, "(n,Xgamma)"},
  {N_XP, "(n,Xp)"},
  {N_XD, "(n,Xd)"},
  {N_XT, "(n,Xt)"},
  {N_X3HE, "(n,X3He)"},
  {N_XA, "(n,Xa)"},
  {HEATING, "heating"},
  {DAMAGE_ENERGY, "damage-energy"},
  {COHERENT, "coherent-scatter"},
  {INCOHERENT, "incoherent-scatter"},
  {PAIR_PROD_ELEC, "pair-production-electron"},
  {PAIR_PROD, "pair-production"},
  {PAIR_PROD_NUC, "pair-production-nuclear"},
  {PHOTOELECTRIC, "photoelectric"},
  {N_PC, "(n,pc)"},
  {N_DC, "(n,dc)"},
  {N_TC, "(n,tc)"},
  {N_3HEC, "(n,3Hec)"},
  {N_AC, "(n,ac)"},
  {N_2NC, "(n,2nc)"},
  {HEATING_LOCAL, "heating-local"},
};

std::string reaction_name(int mt)
{
  if (N_N1 <= mt && mt <= N_N40) {
    return fmt::format("(n,n{})", mt - 50);
  } else if (534 <= mt && mt <= 572) {
    return fmt::format("photoelectric, {} subshell", SUBSHELLS[mt - 534]);
  } else if (N_P0 <= mt && mt < N_PC) {
    return fmt::format("(n,p{})", mt - N_P0);
  } else if (N_D0 <= mt && mt < N_DC) {
    return fmt::format("(n,d{})", mt - N_D0);
  } else if (N_T0 <= mt && mt < N_TC) {
    return fmt::format("(n,t{})", mt - N_T0);
  } else if (N_3HE0 <= mt && mt < N_3HEC) {
    return fmt::format("(n,3He{})", mt - N_3HE0);
  } else if (N_A0 <= mt && mt < N_AC) {
    return fmt::format("(n,a{})", mt - N_A0);
  } else if (N_2N0 <= mt && mt < N_2NC) {
    return fmt::format("(n,2n{})", mt - N_2N0);
  } else {
    auto it = REACTION_NAME_MAP.find(mt);
    if (it != REACTION_NAME_MAP.end()) {
      return it->second;
    } else {
      return fmt::format("MT={}", mt);
    }
  }
}

} // namespace openmc
