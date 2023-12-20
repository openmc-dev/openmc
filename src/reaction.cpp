#include "openmc/reaction.h"

#include <string>
#include <unordered_map>
#include <utility> // for move

#include <fmt/core.h>

#include "openmc/constants.h"
#include "openmc/endf.h"
#include "openmc/hdf5_interface.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/secondary_uncorrelated.h"

namespace openmc {

//==============================================================================
// Reaction implementation
//==============================================================================

Reaction::Reaction(hid_t group, const vector<int>& temperatures)
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
}

double Reaction::xs(
  gsl::index i_temp, gsl::index i_grid, double interp_factor) const
{
  // If energy is below threshold, return 0. Otherwise interpolate between
  // nearest grid points
  const auto& x = xs_[i_temp];
  return (i_grid < x.threshold)
           ? 0.0
           : (1.0 - interp_factor) * x.value[i_grid - x.threshold] +
               interp_factor * x.value[i_grid - x.threshold + 1];
}

double Reaction::xs(const NuclideMicroXS& micro) const
{
  return this->xs(micro.index_temp, micro.index_grid, micro.interp_factor);
}

double Reaction::collapse_rate(gsl::index i_temp,
  gsl::span<const double> energy, gsl::span<const double> flux,
  const vector<double>& grid) const
{
  // Find index corresponding to first energy
  const auto& xs = xs_[i_temp].value;
  int i_low = lower_bound_index(grid.cbegin(), grid.cend(), energy.front());

  // Check for threshold and adjust starting point if necessary
  int j_start = 0;
  int i_threshold = xs_[i_temp].threshold;
  if (i_low < i_threshold) {
    i_low = i_threshold;
    while (energy[j_start + 1] < grid[i_low]) {
      ++j_start;
      if (j_start + 1 == energy.size())
        return 0.0;
    }
  }

  double xs_flux_sum = 0.0;

  for (int j = j_start; j < flux.size(); ++j) {
    double E_group_low = energy[j];
    double E_group_high = energy[j + 1];
    double flux_per_eV = flux[j] / (E_group_high - E_group_low);

    // Determine energy grid index corresponding to group high
    int i_high = i_low;
    while (grid[i_high + 1] < E_group_high && i_high + 1 < grid.size() - 1)
      ++i_high;

    // Loop over energy grid points within [E_group_low, E_group_high]
    for (; i_low <= i_high; ++i_low) {
      // Determine bounding grid energies and cross sections
      double E_l = grid[i_low];
      double E_r = grid[i_low + 1];
      if (E_l == E_r)
        continue;

      double xs_l = xs[i_low - i_threshold];
      double xs_r = xs[i_low + 1 - i_threshold];

      // Determine actual energies
      double E_low = std::max(E_group_low, E_l);
      double E_high = std::min(E_group_high, E_r);

      // Determine average cross section across segment
      double m = (xs_r - xs_l) / (E_r - E_l);
      double xs_low = xs_l + m * (E_low - E_l);
      double xs_high = xs_l + m * (E_high - E_l);
      double xs_avg = 0.5 * (xs_low + xs_high);

      // Add contribution from segment
      double dE = (E_high - E_low);
      xs_flux_sum += flux_per_eV * xs_avg * dE;
    }

    i_low = i_high;

    // Check for end of energy grid
    if (i_low + 1 == grid.size())
      break;
  }

  return xs_flux_sum;
}

//==============================================================================
// Non-member functions
//==============================================================================

std::unordered_map<int, std::string> REACTION_NAME_MAP {
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
  {SCORE_PULSE_HEIGHT, "pulse-height"},
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
  {N_5N, "(n,5n)"},
  {N_6N, "(n,6n)"},
  {N_2NT, "(n,2nt)"},
  {N_TA, "(n,ta)"},
  {N_4NP, "(n,4np)"},
  {N_3ND, "(n,3nd)"},
  {N_NDA, "(n,nda)"},
  {N_2NPA, "(n,2npa)"},
  {N_7N, "(n,7n)"},
  {N_8N, "(n,8n)"},
  {N_5NP, "(n,5np)"},
  {N_6NP, "(n,6np)"},
  {N_7NP, "(n,7np)"},
  {N_4NA, "(n,4na)"},
  {N_5NA, "(n,5na)"},
  {N_6NA, "(n,6na)"},
  {N_7NA, "(n,7na)"},
  {N_4ND, "(n,4nd)"},
  {N_5ND, "(n,5nd)"},
  {N_6ND, "(n,6nd)"},
  {N_3NT, "(n,3nt)"},
  {N_4NT, "(n,4nt)"},
  {N_5NT, "(n,5nt)"},
  {N_6NT, "(n,6nt)"},
  {N_2N3HE, "(n,2n3He)"},
  {N_3N3HE, "(n,3n3He)"},
  {N_4N3HE, "(n,4n3He)"},
  {N_3N2P, "(n,3n2p)"},
  {N_3N2A, "(n,3n2a)"},
  {N_3NPA, "(n,3npa)"},
  {N_DT, "(n,dt)"},
  {N_NPD, "(n,npd)"},
  {N_NPT, "(n,npt)"},
  {N_NDT, "(n,ndt)"},
  {N_NP3HE, "(n,np3He)"},
  {N_ND3HE, "(n,nd3He)"},
  {N_NT3HE, "(n,nt3He)"},
  {N_NTA, "(n,nta)"},
  {N_2N2P, "(n,2n2p)"},
  {N_P3HE, "(n,p3He)"},
  {N_D3HE, "(n,d3He)"},
  {N_3HEA, "(n,3Hea)"},
  {N_4N2P, "(n,4n2p)"},
  {N_4N2A, "(n,4n2a)"},
  {N_4NPA, "(n,4npa)"},
  {N_3P, "(n,3p)"},
  {N_N3P, "(n,n3p)"},
  {N_3N2PA, "(n,3n2pa)"},
  {N_5N2P, "(n,5n2p)"},
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

std::unordered_map<std::string, int> REACTION_TYPE_MAP;

void initialize_maps()
{
  // Add level reactions to name map
  for (int level = 0; level <= 48; ++level) {
    if (level >= 1 && level <= 40) {
      REACTION_NAME_MAP[50 + level] = fmt::format("(n,n{})", level);
    }
    REACTION_NAME_MAP[600 + level] = fmt::format("(n,p{})", level);
    REACTION_NAME_MAP[650 + level] = fmt::format("(n,d{})", level);
    REACTION_NAME_MAP[700 + level] = fmt::format("(n,t{})", level);
    REACTION_NAME_MAP[750 + level] = fmt::format("(n,3He{})", level);
    REACTION_NAME_MAP[800 + level] = fmt::format("(n,a{})", level);
    if (level <= 15) {
      REACTION_NAME_MAP[875 + level] = fmt::format("(n,2n{})", level);
    }
  }

  // Create photoelectric subshells
  for (int mt = 534; mt <= 572; ++mt) {
    REACTION_NAME_MAP[mt] =
      fmt::format("photoelectric, {} subshell", SUBSHELLS[mt - 534]);
  }

  // Invert name map to create type map
  for (const auto& kv : REACTION_NAME_MAP) {
    REACTION_TYPE_MAP[kv.second] = kv.first;
  }
}

std::string reaction_name(int mt)
{
  // Initialize remainder of name map and all of type map
  if (REACTION_TYPE_MAP.empty())
    initialize_maps();

  // Get reaction name from map
  auto it = REACTION_NAME_MAP.find(mt);
  if (it != REACTION_NAME_MAP.end()) {
    return it->second;
  } else {
    return fmt::format("MT={}", mt);
  }
}

int reaction_type(std::string name)
{
  // Initialize remainder of name map and all of type map
  if (REACTION_TYPE_MAP.empty())
    initialize_maps();

  // (n,total) exists in REACTION_TYPE_MAP for MT=1, but we need this to return
  // the special SCORE_TOTAL score
  if (name == "(n,total)")
    return SCORE_TOTAL;

  // Check if type map has an entry for this reaction name
  auto it = REACTION_TYPE_MAP.find(name);
  if (it != REACTION_TYPE_MAP.end()) {
    return it->second;
  }

  // Alternate names for several reactions
  if (name == "elastic") {
    return ELASTIC;
  } else if (name == "n2n") {
    return N_2N;
  } else if (name == "n3n") {
    return N_3N;
  } else if (name == "n4n") {
    return N_4N;
  } else if (name == "H1-production") {
    return N_XP;
  } else if (name == "H2-production") {
    return N_XD;
  } else if (name == "H3-production") {
    return N_XT;
  } else if (name == "He3-production") {
    return N_X3HE;
  } else if (name == "He4-production") {
    return N_XA;
  }

  // Assume the given string is a reaction MT number.  Make sure it's a natural
  // number then return.
  int MT = 0;
  try {
    MT = std::stoi(name);
  } catch (const std::invalid_argument& ex) {
    throw std::invalid_argument(
      "Invalid tally score \"" + name +
      "\". See the docs "
      "for details: "
      "https://docs.openmc.org/en/stable/usersguide/tallies.html#scores");
  }
  if (MT < 1)
    throw std::invalid_argument(
      "Invalid tally score \"" + name +
      "\". See the docs "
      "for details: "
      "https://docs.openmc.org/en/stable/usersguide/tallies.html#scores");
  return MT;
}

} // namespace openmc
