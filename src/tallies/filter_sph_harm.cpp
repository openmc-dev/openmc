#include "openmc/tallies/filter_sph_harm.h"

#include <utility>  // For pair

#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
SphericalHarmonicsFilter::from_xml(pugi::xml_node node)
{
  order_ = std::stoi(get_node_value(node, "order"));
  n_bins_ = (order_ + 1) * (order_ + 1);

  if (check_for_node(node, "cosine")) {
    auto cos = get_node_value(node, "cosine", true);
    if (cos == "scatter") {
      cosine_ = SphericalHarmonicsCosine::scatter;
    } else if (cos == "particle") {
      cosine_ = SphericalHarmonicsCosine::particle;
    } else {
      std::stringstream err_msg;
      err_msg << "Unrecognized cosine type, \"" << cos
              << "\" in spherical harmonics filter";
      fatal_error(err_msg);
    }
  }
}

void
SphericalHarmonicsFilter::get_all_bins(const Particle* p, int estimator,
                                       FilterMatch& match) const
{
  // Determine cosine term for scatter expansion if necessary
  std::vector<double> wgt(order_ + 1);
  if (cosine_ == SphericalHarmonicsCosine::scatter) {
    calc_pn_c(order_, p->mu_, wgt.data());
  } else {
    for (int i = 0; i < order_ + 1; i++) wgt[i] = 1;
  }

  // Find the Rn,m values
  std::vector<double> rn(n_bins_);
  calc_rn(order_, p->u_last_, rn.data());

  int j = 0;
  for (int n = 0; n < order_ + 1; n++) {
    // Calculate n-th order spherical harmonics for (u,v,w)
    int num_nm = 2*n + 1;

    // Append the matching (bin,weight) for each moment
    for (int i = 0; i < num_nm; i++) {
      match.weights_.push_back(wgt[n] * rn[j]);
      match.bins_.push_back(j);
      ++j;
    }
  }
}

void
SphericalHarmonicsFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "order", order_);
  if (cosine_ == SphericalHarmonicsCosine::scatter) {
    write_dataset(filter_group, "cosine", "scatter");
  } else {
    write_dataset(filter_group, "cosine", "particle");
  }
}

std::string
SphericalHarmonicsFilter::text_label(int bin) const
{
  std::stringstream out;
  for (int n = 0; n < order_ + 1; n++) {
    if (bin < (n + 1) * (n + 1)) {
      int m = (bin - n*n) - n;
      out << "Spherical harmonic expansion, Y" << n << "," << m;
      break;
    }
  }
  return out.str();
}

//==============================================================================
// C-API functions
//==============================================================================

std::pair<int, SphericalHarmonicsFilter*>
check_sphharm_filter(int32_t index)
{
  // Make sure this is a valid index to an allocated filter.
  int err = verify_filter(index);
  if (err) {
    return {err, nullptr};
  }

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<SphericalHarmonicsFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a spherical harmonics filter.");
    err = OPENMC_E_INVALID_TYPE;
  }
  return {err, filt};
}

extern "C" int
openmc_sphharm_filter_get_order(int32_t index, int* order)
{
  // Check the filter.
  auto check_result = check_sphharm_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Output the order.
  *order = filt->order_;
  return 0;
}

extern "C" int
openmc_sphharm_filter_get_cosine(int32_t index, char cosine[])
{
  // Check the filter.
  auto check_result = check_sphharm_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Output the cosine.
  if (filt->cosine_ == SphericalHarmonicsCosine::scatter) {
    strcpy(cosine, "scatter");
  } else {
    strcpy(cosine, "particle");
  }
  return 0;
}

extern "C" int
openmc_sphharm_filter_set_order(int32_t index, int order)
{
  // Check the filter.
  auto check_result = check_sphharm_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Update the filter.
  filt->order_ = order;
  filt->n_bins_ = (order + 1) * (order + 1);
  return 0;
}

extern "C" int
openmc_sphharm_filter_set_cosine(int32_t index, const char cosine[])
{
  // Check the filter.
  auto check_result = check_sphharm_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Update the filter.
  if (strcmp(cosine, "scatter") == 0) {
    filt->cosine_ = SphericalHarmonicsCosine::scatter;
  } else if (strcmp(cosine, "particle") == 0) {
    filt->cosine_ = SphericalHarmonicsCosine::particle;
  } else {
    set_errmsg("Invalid spherical harmonics cosine.");
    return OPENMC_E_INVALID_ARGUMENT;
  }
  return 0;
}

} // namespace openmc
