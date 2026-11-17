#include "openmc/secondary_data.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "openmc/chain.h"
#include "openmc/distribution.h"
#include "openmc/distribution_angle.h"
#include "openmc/distribution_energy.h"
#include "openmc/endf.h"
#include "openmc/error.h"
#include "openmc/memory.h"
#include "openmc/nuclide.h"
#include "openmc/random_dist.h"
#include "openmc/reaction.h"
#include "openmc/reaction_product.h"
#include "openmc/secondary_correlated.h"
#include "openmc/secondary_kalbach.h"
#include "openmc/secondary_nbody.h"
#include "openmc/secondary_uncorrelated.h"
#include "openmc/simulation.h"

namespace openmc {

Nuclide& secondary_nuclide(int index)
{
  if (simulation::in_batch)
    throw std::invalid_argument(
      "Secondary data cannot be inspected during a batch.");
  if (index < 0 || index >= data::nuclides.size())
    throw std::out_of_range("Loaded nuclide index is out of range.");
  return *data::nuclides[index];
}

Reaction& secondary_reaction(int index, int mt)
{
  auto& nuc = secondary_nuclide(index);
  if (mt < 0 || mt >= nuc.reaction_index_.size() ||
      nuc.reaction_index_[mt] == C_NONE)
    throw std::out_of_range("No reaction with requested MT number.");
  return *nuc.reactions_[nuc.reaction_index_[mt]];
}

ReactionProduct& secondary_product(int index, int mt, int product)
{
  auto& rx = secondary_reaction(index, mt);
  if (product < 0 || product >= rx.products_.size())
    throw std::out_of_range("Reaction product index is out of range.");
  return rx.products_[product];
}

AngleEnergy& secondary_distribution(int index, int mt, int product, int dist)
{
  auto& p = secondary_product(index, mt, product);
  if (dist < 0 || dist >= p.distribution_.size())
    throw std::out_of_range("Product distribution index is out of range.");
  if (!p.distribution_[dist])
    throw std::runtime_error("Null native secondary distribution.");
  return *p.distribution_[dist];
}

// Inspection adapters have private access but accept only const native objects.
// No adapter evaluates, normalizes, samples, or changes the native data.
class SecondaryDataReader {
public:
  SecondarySnapshot snapshot;

  size_t node(size_t parent, const std::string& key, int type)
  {
    size_t id = snapshot.nodes.size();
    snapshot.nodes.emplace_back();
    snapshot.nodes[id].name = key;
    snapshot.nodes[id].type = type;
    if (id != 0)
      snapshot.nodes[parent].children.push_back(id);
    return id;
  }

  size_t object(
    size_t parent, const std::string& key, const char* type = nullptr)
  {
    auto id = node(parent, key, OPENMC_SECONDARY_OBJECT);
    if (type)
      text(id, "type", type);
    return id;
  }

  void text(size_t parent, const std::string& key, const std::string& value)
  {
    auto id = node(parent, key, OPENMC_SECONDARY_STRING);
    snapshot.nodes[id].text = value;
  }

  void integer(size_t parent, const std::string& key, int64_t value,
    int type = OPENMC_SECONDARY_INT)
  {
    auto id = node(parent, key, type);
    snapshot.nodes[id].integers = {value};
  }

  void real(size_t parent, const std::string& key, double value)
  {
    auto id = node(parent, key, OPENMC_SECONDARY_DOUBLE);
    snapshot.nodes[id].reals = {value};
  }

  template<class T>
  void reals(size_t parent, const std::string& key, const T& values)
  {
    auto id = node(parent, key, OPENMC_SECONDARY_DOUBLE_ARRAY);
    snapshot.nodes[id].reals.assign(values.begin(), values.end());
  }

  template<class T>
  void integers(
    size_t parent, const std::string& key, const T& values, int offset = 0)
  {
    auto id = node(parent, key, OPENMC_SECONDARY_INT_ARRAY);
    for (auto value : values)
      snapshot.nodes[id].integers.push_back(
        static_cast<int64_t>(value) + offset);
  }

  void function(size_t parent, const std::string& key, const Function1D& f)
  {
    auto id = object(parent, key);
    if (auto p = dynamic_cast<const Tabulated1D*>(&f)) {
      text(id, "type", "tabulated-1d");
      reals(id, "x", p->x_);
      reals(id, "y", p->y_);
      // Public breakpoint convention matches ENDF/HDF5 (one-based endpoints).
      integers(id, "breakpoints", p->nbt_, 1);
      integers(id, "interpolation", p->int_);
    } else if (auto p = dynamic_cast<const Polynomial*>(&f)) {
      text(id, "type", "polynomial");
      reals(id, "coefficients", p->coef_);
    } else if (auto p = dynamic_cast<const CoherentElasticXS*>(&f)) {
      text(id, "type", "coherent-elastic-xs");
      reals(id, "bragg_edges", p->bragg_edges_);
      reals(id, "factors", p->factors_);
    } else if (auto p = dynamic_cast<const IncoherentElasticXS*>(&f)) {
      text(id, "type", "incoherent-elastic-xs");
      real(id, "bound_xs", p->bound_xs_);
      real(id, "debye_waller", p->debye_waller_);
    } else if (auto p = dynamic_cast<const Sum1D*>(&f)) {
      text(id, "type", "sum");
      auto list = node(id, "functions", OPENMC_SECONDARY_LIST);
      for (const auto& child : p->functions_)
        function(list, "", *child);
    } else {
      throw std::runtime_error(
        "Unsupported native Function1D inspection adapter.");
    }
  }

  void tabular(size_t parent, const std::string& key, const Tabular& d)
  {
    auto id = object(parent, key, "tabular");
    reals(id, "x", d.x_);
    reals(id, "pdf", d.p_);
    reals(id, "cdf", d.c_);
    integer(id, "interpolation", static_cast<int>(d.interp_));
    real(id, "integral", d.integral_);
  }

  void discrete_index(size_t parent, const DiscreteIndex& d)
  {
    reals(parent, "alias_probability", d.prob());
    integers(parent, "alias", d.alias());
    real(parent, "alias_integral", d.integral());
  }

  // Chain spectra use generic Distribution objects, not EnergyDistribution.
  // Preserve native alias tables, mixture weights, and bias information too.
  void generic(size_t parent, const std::string& key, const Distribution& d)
  {
    if (auto p = dynamic_cast<const Tabular*>(&d)) {
      tabular(parent, key, *p);
      auto id = snapshot.nodes[parent].children.back();
      if (d.bias())
        generic(id, "bias", *d.bias());
      return;
    }
    auto id = object(parent, key);
    real(id, "integral", d.integral());
    if (auto p = dynamic_cast<const Discrete*>(&d)) {
      text(id, "type", "discrete");
      reals(id, "x", p->x_);
      reals(id, "importance_weights", p->weight_);
      discrete_index(id, p->di_);
    } else if (auto p = dynamic_cast<const Uniform*>(&d)) {
      text(id, "type", "uniform");
      real(id, "a", p->a_);
      real(id, "b", p->b_);
    } else if (auto p = dynamic_cast<const PowerLaw*>(&d)) {
      text(id, "type", "power-law");
      real(id, "offset", p->offset_);
      real(id, "span", p->span_);
      real(id, "inverse_exponent", p->ninv_);
    } else if (auto p = dynamic_cast<const Maxwell*>(&d)) {
      text(id, "type", "maxwell");
      real(id, "theta", p->theta_);
    } else if (auto p = dynamic_cast<const Watt*>(&d)) {
      text(id, "type", "watt");
      real(id, "a", p->a_);
      real(id, "b", p->b_);
    } else if (auto p = dynamic_cast<const Normal*>(&d)) {
      text(id, "type", "normal");
      real(id, "mean", p->mean_value_);
      real(id, "std_dev", p->std_dev_);
      real(id, "lower", p->lower_);
      real(id, "upper", p->upper_);
      real(id, "normalization", p->norm_factor_);
      integer(id, "truncated", p->is_truncated_, OPENMC_SECONDARY_BOOL);
    } else if (auto p = dynamic_cast<const Equiprobable*>(&d)) {
      text(id, "type", "equiprobable");
      reals(id, "x", p->x_);
    } else if (auto p = dynamic_cast<const Mixture*>(&d)) {
      text(id, "type", "mixture");
      discrete_index(id, p->di_);
      reals(id, "importance_weights", p->weight_);
      auto list = node(id, "components", OPENMC_SECONDARY_LIST);
      for (const auto& child : p->distribution_)
        generic(list, "", *child);
    } else if (auto p = dynamic_cast<const DecaySpectrum*>(&d)) {
      text(id, "type", "decay-spectrum");
      discrete_index(id, p->di_);
      integers(id, "chain_nuclide_indices", p->nuclide_indices_);
      auto list = node(id, "components", OPENMC_SECONDARY_LIST);
      for (auto index : p->nuclide_indices_)
        generic(list, "", *data::chain_nuclides.at(index)->photon_energy());
    } else {
      throw std::runtime_error(
        "Unsupported chain spectrum inspection adapter.");
    }
    if (d.bias())
      generic(id, "bias", *d.bias());
  }

  template<class T>
  void incident(size_t id, const T& d)
  {
    reals(id, "incident_energy", d.energy_);
    integers(id, "breakpoints", d.breakpoints_);
    integers(id, "interpolation", d.interpolation_);
  }

  template<class T>
  size_t energy_table(size_t list, const T& table)
  {
    auto id = object(list, "");
    reals(id, "energy", table.e_out);
    reals(id, "pdf", table.p);
    reals(id, "cdf", table.c);
    integer(id, "interpolation", static_cast<int>(table.interpolation));
    integer(id, "n_discrete", table.n_discrete);
    return id;
  }

  static const char* energy_type(const EnergyDistribution* d)
  {
    if (!d)
      return nullptr;
    if (dynamic_cast<const ContinuousTabular*>(d))
      return "continuous-tabular";
    if (dynamic_cast<const DiscretePhoton*>(d))
      return "discrete-photon";
    if (dynamic_cast<const LevelInelastic*>(d))
      return "level-inelastic";
    if (dynamic_cast<const MaxwellEnergy*>(d))
      return "maxwell";
    if (dynamic_cast<const Evaporation*>(d))
      return "evaporation";
    if (dynamic_cast<const WattEnergy*>(d))
      return "watt";
    return "unsupported";
  }

  static const char* angle_energy_type(const AngleEnergy& d)
  {
    if (dynamic_cast<const UncorrelatedAngleEnergy*>(&d))
      return "uncorrelated";
    if (dynamic_cast<const CorrelatedAngleEnergy*>(&d))
      return "correlated";
    if (dynamic_cast<const KalbachMann*>(&d))
      return "kalbach-mann";
    if (dynamic_cast<const NBodyPhaseSpace*>(&d))
      return "nbody";
    if (dynamic_cast<const DecayPhotonAngleEnergy*>(&d))
      return "decay-photon";
    return "unsupported";
  }

  void energy(size_t parent, const EnergyDistribution* d)
  {
    if (!d) {
      node(parent, "energy", OPENMC_SECONDARY_NULL);
      return;
    }
    auto id = object(parent, "energy", energy_type(d));
    if (auto p = dynamic_cast<const ContinuousTabular*>(d)) {
      incident(id, *p);
      text(id, "sampling_interpolation",
        p->n_region_ == 1 && p->interpolation_[0] == Interpolation::histogram
          ? "histogram"
          : "linear-stochastic-unit-base");
      auto list = node(id, "tables", OPENMC_SECONDARY_LIST);
      for (const auto& table : p->distribution_)
        energy_table(list, table);
    } else if (auto p = dynamic_cast<const DiscretePhoton*>(d)) {
      integer(id, "primary_flag", p->primary_flag_);
      real(id, "energy", p->energy_);
      real(id, "atomic_weight_ratio", p->A_);
    } else if (auto p = dynamic_cast<const LevelInelastic*>(d)) {
      real(id, "threshold", p->threshold_);
      real(id, "mass_ratio", p->mass_ratio_);
    } else if (auto p = dynamic_cast<const MaxwellEnergy*>(d)) {
      function(id, "theta", p->theta_);
      real(id, "u", p->u_);
    } else if (auto p = dynamic_cast<const Evaporation*>(d)) {
      function(id, "theta", p->theta_);
      real(id, "u", p->u_);
    } else if (auto p = dynamic_cast<const WattEnergy*>(d)) {
      function(id, "a", p->a_);
      function(id, "b", p->b_);
      real(id, "u", p->u_);
    } else {
      throw std::runtime_error(
        "Unsupported EnergyDistribution inspection adapter.");
    }
  }

  void angle(size_t parent, const AngleDistribution& d)
  {
    auto id =
      object(parent, "angle", d.empty() ? "isotropic" : "tabular-angle");
    if (d.empty())
      return;
    reals(id, "incident_energy", d.energy_);
    text(id, "sampling_interpolation", "linear-stochastic");
    auto list = node(id, "tables", OPENMC_SECONDARY_LIST);
    for (const auto& table : d.distribution_)
      tabular(list, "", *table);
  }

  void law(size_t parent, const std::string& key, const AngleEnergy& d)
  {
    auto id = object(parent, key, angle_energy_type(d));
    if (auto p = dynamic_cast<const UncorrelatedAngleEnergy*>(&d)) {
      angle(id, p->angle_);
      energy(id, p->energy_.get());
    } else if (auto p = dynamic_cast<const CorrelatedAngleEnergy*>(&d)) {
      incident(id, *p);
      text(id, "sampling_interpolation", "linear-stochastic-unit-base");
      auto list = node(id, "tables", OPENMC_SECONDARY_LIST);
      for (const auto& table : p->distribution_) {
        auto t = energy_table(list, table);
        auto angles = node(t, "angle", OPENMC_SECONDARY_LIST);
        for (const auto& mu : table.angle)
          tabular(angles, "", *mu);
      }
    } else if (auto p = dynamic_cast<const KalbachMann*>(&d)) {
      incident(id, *p);
      text(id, "sampling_interpolation", "linear-stochastic-unit-base");
      auto list = node(id, "tables", OPENMC_SECONDARY_LIST);
      for (const auto& table : p->distribution_) {
        auto t = energy_table(list, table);
        reals(t, "r", table.r);
        reals(t, "a", table.a);
      }
    } else if (auto p = dynamic_cast<const NBodyPhaseSpace*>(&d)) {
      object(id, "angle", "isotropic");
      integer(id, "n_particles", p->n_bodies_);
      real(id, "total_mass", p->mass_ratio_);
      real(id, "atomic_weight_ratio", p->A_);
      real(id, "q_value", p->Q_);
    } else if (auto p = dynamic_cast<const DecayPhotonAngleEnergy*>(&d)) {
      object(id, "angle", "isotropic");
      text(id, "ownership", "borrowed-from-chain-nuclide");
      if (!p->photon_energy_)
        throw std::runtime_error("Missing chain photon distribution.");
      generic(id, "spectrum", *p->photon_energy_);
    } else {
      throw std::runtime_error("Unsupported AngleEnergy inspection adapter.");
    }
  }

  void product(size_t root, const Nuclide& nuc, const Reaction& rx, int ip,
    int selected_dist, bool include_data)
  {
    const auto& p = rx.products_[ip];
    auto id = object(root, "");
    integer(id, "schema_version", 1);
    text(id, "nuclide", nuc.name_);
    integer(id, "nuclide_index", nuc.index_);
    integer(id, "mt", rx.mt_);
    integer(id, "product_index", ip);
    text(id, "particle", p.particle_.str());
    const char* mode =
      p.emission_mode_ == ReactionProduct::EmissionMode::prompt    ? "prompt"
      : p.emission_mode_ == ReactionProduct::EmissionMode::delayed ? "delayed"
                                                                   : "total";
    text(id, "emission_mode", mode);
    text(id, "reference_frame",
      rx.scatter_in_cm_ ? "center-of-mass" : "laboratory");
    integer(id, "parent_nuclide_index", p.parent_nuclide_);
    if (p.yield_)
      function(id, "yield", *p.yield_);
    else
      node(id, "yield", OPENMC_SECONDARY_NULL);
    auto apps = node(id, "applicability", OPENMC_SECONDARY_LIST);
    for (const auto& app : p.applicability_)
      function(apps, "", app);
    auto list = node(id, "distributions", OPENMC_SECONDARY_LIST);
    for (int i = 0; i < p.distribution_.size(); ++i) {
      if (selected_dist != -1 && selected_dist != i)
        continue;
      const auto& d = *p.distribution_[i];
      auto entry = object(list, "");
      integer(entry, "distribution_index", i);
      text(entry, "angle_energy_type", angle_energy_type(d));
      auto uncorrelated = dynamic_cast<const UncorrelatedAngleEnergy*>(&d);
      const auto* e = uncorrelated ? uncorrelated->energy() : nullptr;
      if (e)
        text(entry, "energy_type", energy_type(e));
      else
        node(entry, "energy_type", OPENMC_SECONDARY_NULL);
      bool known = std::string(angle_energy_type(d)) != "unsupported" &&
                   (!e || std::string(energy_type(e)) != "unsupported");
      integer(entry, "read_supported", known, OPENMC_SECONDARY_BOOL);
      integer(entry, "mutation_supported", false, OPENMC_SECONDARY_BOOL);
      bool angular = uncorrelated && !uncorrelated->angle().empty();
      if (angular) {
        for (const auto& table : uncorrelated->angle_.distribution_)
          angular = angular && (table->interp_ == Interpolation::histogram ||
                                 table->interp_ == Interpolation::lin_lin);
      }
      integer(
        entry, "angular_mutation_supported", angular, OPENMC_SECONDARY_BOOL);
      text(entry, "angular_mutation_reason",
        angular ? "Existing histogram/linear-linear tabular angular PDFs on "
                  "unchanged grids."
                : "Isotropic representation changes and conditional angular "
                  "mutations are not implemented.");
      bool ct = dynamic_cast<const ContinuousTabular*>(e) != nullptr;
      integer(entry, "energy_mutation_supported", ct, OPENMC_SECONDARY_BOOL);
      if (ct) {
        text(entry, "mutation_reason",
          "Whole-distribution mutation is not implemented; "
          "energy-only editing is available through set_emission_spectrum.");
        text(entry, "energy_mutation_api", "set_emission_spectrum");
      } else if (dynamic_cast<const DecayPhotonAngleEnergy*>(&d)) {
        text(entry, "mutation_reason",
          "Read-only: spectrum is borrowed from a chain "
          "nuclide and may be shared by other products and consumers.");
      } else {
        text(entry, "mutation_reason",
          "Runtime mutation adapter is not implemented.");
      }
      if (!known)
        text(entry, "read_reason",
          "No native inspection adapter for this representation.");
      if (include_data && known)
        law(entry, "data", d);
    }
  }
};

SecondarySnapshot inspect_secondary(const AngleEnergy& distribution)
{
  SecondaryDataReader reader;
  reader.law(0, "", distribution);
  return std::move(reader.snapshot);
}

// Mutation adapter owns candidates until every table has passed validation.
class SecondaryDataWriter {
public:
  static unique_ptr<Tabular> clone_table(const Tabular& original)
  {
    if (original.bias())
      throw std::invalid_argument("Biased angular tables cannot be copied.");
    auto table = make_unique<Tabular>(original.x_.data(), original.p_.data(),
      original.x_.size(), original.interp_, original.c_.data());
    table->p_ = original.p_;
    table->c_ = original.c_;
    table->integral_ = original.integral_;
    return table;
  }

  static AngleDistribution clone(const AngleDistribution& angle)
  {
    AngleDistribution result;
    result.energy_ = angle.energy_;
    for (const auto& original : angle.distribution_) {
      if (original->bias())
        throw std::invalid_argument(
          "Biased angular tables cannot be backed up.");
      auto table =
        make_unique<Tabular>(original->x_.data(), original->p_.data(),
          original->x_.size(), original->interp_, original->c_.data());
      // Bypass reconstruction/renormalization: preserve loaded values bitwise.
      table->p_ = original->p_;
      table->c_ = original->c_;
      table->integral_ = original->integral_;
      result.distribution_.push_back(std::move(table));
    }
    return result;
  }

  static void set_angle(AngleDistribution& angle, const double* incident, int n,
    const int* offsets, const int* interpolation, const double* mu,
    const double* pdf, int count)
  {
    if (!incident || !offsets || !interpolation || !mu || !pdf || n <= 0 ||
        count <= 0 || n != angle.energy_.size() || offsets[0] != 0 ||
        offsets[n] != count)
      throw std::invalid_argument(
        "Invalid angular arrays or native table count.");
    vector<unique_ptr<Tabular>> tables;
    for (int i = 0; i < n; ++i) {
      const auto& original = *angle.distribution_[i];
      int lo = offsets[i], hi = offsets[i + 1];
      if (incident[i] != angle.energy_[i] || lo < 0 || hi > count || hi <= lo ||
          hi - lo != original.x_.size() ||
          interpolation[i] != static_cast<int>(original.interp_) ||
          (interpolation[i] != 1 && interpolation[i] != 2))
        throw std::invalid_argument(
          "Angular grids and interpolation must match native data.");
      for (int j = lo; j < hi; ++j) {
        if (!std::isfinite(mu[j]) || mu[j] < -1.0 || mu[j] > 1.0 ||
            mu[j] != original.x_[j - lo] || !std::isfinite(pdf[j]) ||
            pdf[j] < 0.0)
          throw std::invalid_argument(
            "Invalid angular coordinate or probability.");
      }
      if (std::equal(pdf + lo, pdf + hi, original.p_.begin())) {
        // Do not reconstruct out-of-window PDFs/CDFs during a partial edit.
        tables.push_back(clone_table(original));
        continue;
      }
      double mass = 0.0;
      for (int j = lo + 1; j < hi; ++j) {
        double width = mu[j] - mu[j - 1];
        if (width < 0.0)
          throw std::invalid_argument("Angular coordinates must be ordered.");
        mass += width * (interpolation[i] == 1 ? pdf[j - 1]
                                               : .5 * (pdf[j - 1] + pdf[j]));
      }
      if (!std::isfinite(mass) || mass <= 0.0)
        throw std::invalid_argument("Invalid angular normalization.");
      auto table =
        make_unique<Tabular>(mu + lo, pdf + lo, hi - lo, original.interp_);
      if (!std::isfinite(table->integral_) || table->integral_ <= 0.0)
        throw std::invalid_argument("Invalid angular normalization.");
      for (auto p : table->p_)
        if (!std::isfinite(p))
          throw std::invalid_argument(
            "Invalid normalized angular probability.");
      tables.push_back(std::move(table));
    }
    angle.distribution_.swap(tables);
  }
};

void replace_secondary_angle(AngleDistribution& angle,
  const vector<double>& incident, const vector<int>& offsets,
  const vector<int>& interpolation, const vector<double>& mu,
  const vector<double>& pdf)
{
  if (offsets.size() != incident.size() + 1 ||
      interpolation.size() != incident.size() || pdf.size() != mu.size())
    throw std::invalid_argument("Invalid angular array lengths.");
  SecondaryDataWriter::set_angle(angle, incident.data(), incident.size(),
    offsets.data(), interpolation.data(), mu.data(), pdf.data(), mu.size());
}

namespace {
struct NuclearDataBackup {
  int index, mt, product, distribution, kind;
  vector<Reaction::TemperatureXS> xs;
  vector<tensor::Tensor<double>> derived;
  unique_ptr<ContinuousTabular> energy;
  AngleDistribution angle;
};
std::unordered_map<uint64_t, unique_ptr<NuclearDataBackup>> backups;
uint64_t next_backup = 1;
} // namespace

void clear_nuclear_data_backups()
{
  backups.clear();
}

extern "C" int openmc_nuclear_data_reaction_info(
  int index, int mt, double* q_value, bool* redundant, bool* center_of_mass)
{
  try {
    if (!q_value || !redundant || !center_of_mass)
      throw std::invalid_argument("Null reaction metadata output.");
    const auto& reaction = secondary_reaction(index, mt);
    *q_value = reaction.q_value_;
    *redundant = reaction.redundant_;
    *center_of_mass = reaction.scatter_in_cm_;
    return 0;
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_INVALID_ARGUMENT;
  }
}

extern "C" int openmc_nuclear_data_backup(
  int index, int mt, int product, int distribution, int kind, uint64_t* handle)
{
  try {
    if (!handle)
      throw std::invalid_argument("Null backup handle.");
    *handle = 0;
    auto& nuc = secondary_nuclide(index);
    auto& rx = secondary_reaction(index, mt);
    auto b = make_unique<NuclearDataBackup>();
    b->index = index;
    b->mt = mt;
    b->product = product;
    b->distribution = distribution;
    b->kind = kind;
    if (kind == 0) {
      // Smooth-array edits are not a consistent URR/multipole/fission model.
      if (rx.redundant_ || nuc.urr_present_ || nuc.multipole_ ||
          nuc.resonant_ || nuc.fissionable_)
        throw std::invalid_argument(
          "Unified XS mutation excludes redundant reactions, URR, multipole, "
          "resonance and fission nuclides.");
      b->xs = rx.xs_;
      b->derived = nuc.xs_;
    } else {
      auto* law = dynamic_cast<UncorrelatedAngleEnergy*>(
        &secondary_distribution(index, mt, product, distribution));
      if (!law)
        throw std::invalid_argument(
          "Conditional secondary laws are read-only.");
      if (kind == 1) {
        auto* energy = dynamic_cast<ContinuousTabular*>(law->energy());
        if (!energy)
          throw std::invalid_argument(
            "Only ContinuousTabular energy rollback supported.");
        b->energy = make_unique<ContinuousTabular>(*energy);
      } else if (kind == 2 && !law->angle().empty()) {
        b->angle = SecondaryDataWriter::clone(law->angle());
      } else {
        throw std::invalid_argument(
          "Unsupported mutation backup kind/representation.");
      }
    }
    auto id = next_backup++;
    backups.emplace(id, std::move(b));
    *handle = id;
    return 0;
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_INVALID_ARGUMENT;
  }
}

extern "C" int openmc_nuclear_data_restore(uint64_t handle)
{
  try {
    auto it = backups.find(handle);
    if (it == backups.end())
      throw std::invalid_argument("Invalid or expired nuclear-data backup.");
    auto& b = *it->second;
    auto& nuc = secondary_nuclide(b.index);
    auto& rx = secondary_reaction(b.index, b.mt);
    // Moves/swaps restore owned native arrays without numerical processing.
    if (b.kind == 0) {
      rx.xs_.swap(b.xs);
      nuc.xs_.swap(b.derived);
    } else {
      auto& law = dynamic_cast<UncorrelatedAngleEnergy&>(
        secondary_distribution(b.index, b.mt, b.product, b.distribution));
      if (b.kind == 1)
        *dynamic_cast<ContinuousTabular*>(law.energy()) = std::move(*b.energy);
      else
        law.angle() = std::move(b.angle);
    }
    backups.erase(it);
    return 0;
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_INVALID_ARGUMENT;
  }
}

extern "C" int openmc_nuclear_data_validate_xs(int index)
{
  try {
    const auto& nuc = secondary_nuclide(index);
    nuc.validate_derived_xs();
    return 0;
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_INVALID_ARGUMENT;
  }
}

extern "C" void openmc_nuclear_data_backup_free(uint64_t handle)
{
  backups.erase(handle);
}

extern "C" int openmc_secondary_set_angle(int index, int mt, int product,
  int distribution, const double* incident, int n, const int* offsets,
  const int* interpolation, const double* mu, const double* pdf, int count)
{
  try {
    if (!simulation::initialized) {
      set_errmsg(
        "Simulation must be initialized before modifying angular data.");
      return OPENMC_E_ALLOCATE;
    }
    auto& law = secondary_distribution(index, mt, product, distribution);
    auto* independent = dynamic_cast<UncorrelatedAngleEnergy*>(&law);
    if (!independent || independent->angle().empty())
      throw std::invalid_argument(
        "Only existing uncorrelated tabular angular laws are mutable.");
    SecondaryDataWriter::set_angle(independent->angle(), incident, n, offsets,
      interpolation, mu, pdf, count);
    return 0;
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_INVALID_ARGUMENT;
  }
}

extern "C" int openmc_secondary_sample_energy(int index, int mt, int product,
  int distribution, double energy, uint64_t* seed, size_t n, double* outgoing)
{
  try {
    if (!seed || (n && !outgoing) || !std::isfinite(energy) || energy < 0.0)
      throw std::invalid_argument("Invalid energy sampling argument.");
    auto& law = secondary_distribution(index, mt, product, distribution);
    auto* independent = dynamic_cast<UncorrelatedAngleEnergy*>(&law);
    if (!independent || !independent->energy())
      throw std::invalid_argument(
        "Requested law has no independent native energy sampler.");
    for (size_t i = 0; i < n; ++i)
      outgoing[i] = independent->energy()->sample(energy, seed);
    return 0;
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_INVALID_ARGUMENT;
  }
}

extern "C" int openmc_secondary_sample_angle(int index, int mt, int product,
  int distribution, double energy, uint64_t* seed, size_t n, double* mu)
{
  try {
    if (!seed || (n && !mu) || !std::isfinite(energy) || energy < 0.0)
      throw std::invalid_argument("Invalid angular sampling argument.");
    auto& law = secondary_distribution(index, mt, product, distribution);
    auto* independent = dynamic_cast<UncorrelatedAngleEnergy*>(&law);
    bool isotropic = dynamic_cast<NBodyPhaseSpace*>(&law) ||
                     dynamic_cast<DecayPhotonAngleEnergy*>(&law);
    if (!independent && !isotropic)
      throw std::invalid_argument(
        "Angular law is conditional on outgoing energy; inspect the full "
        "emission distribution.");
    // This is the same AngleDistribution::sample used by elastic transport.
    // The seed is local to the caller; transport RNG streams are untouched.
    for (size_t i = 0; i < n; ++i) {
      mu[i] = independent && !independent->angle().empty()
                ? independent->angle().sample(energy, seed)
                : uniform_distribution(-1.0, 1.0, seed);
    }
    return 0;
  } catch (const std::out_of_range& e) {
    set_errmsg(e.what());
    return OPENMC_E_OUT_OF_BOUNDS;
  } catch (const std::invalid_argument& e) {
    set_errmsg(e.what());
    return OPENMC_E_INVALID_ARGUMENT;
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_DATA;
  }
}

extern "C" int openmc_secondary_snapshot_create(int index, int mt, int product,
  int distribution, bool include_data, void** result)
{
  if (!result) {
    set_errmsg("Null secondary snapshot output pointer.");
    return OPENMC_E_INVALID_ARGUMENT;
  }
  *result = nullptr;
  try {
    if (mt < -1 || product < -1 || distribution < -1 ||
        (mt == -1 && product != -1) || (product == -1 && distribution != -1))
      throw std::invalid_argument("Invalid secondary distribution selection.");
    const auto& nuc = secondary_nuclide(index);
    if (mt != -1)
      secondary_reaction(index, mt);
    if (product != -1)
      secondary_product(index, mt, product);
    if (distribution != -1)
      secondary_distribution(index, mt, product, distribution);
    SecondaryDataReader reader;
    auto root = reader.node(0, "", OPENMC_SECONDARY_LIST);
    for (const auto& rx : nuc.reactions_) {
      if (mt != -1 && rx->mt_ != mt)
        continue;
      for (int ip = 0; ip < rx->products_.size(); ++ip) {
        if (product != -1 && product != ip)
          continue;
        reader.product(root, nuc, *rx, ip, distribution, include_data);
      }
    }
    *result = new SecondarySnapshot(std::move(reader.snapshot));
    return 0;
  } catch (const std::out_of_range& e) {
    set_errmsg(e.what());
    return OPENMC_E_OUT_OF_BOUNDS;
  } catch (const std::invalid_argument& e) {
    set_errmsg(e.what());
    return OPENMC_E_INVALID_ARGUMENT;
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_DATA;
  }
}

extern "C" int openmc_secondary_snapshot_node(
  const void* snapshot, size_t index, openmc_secondary_node* result)
{
  if (!snapshot || !result) {
    set_errmsg("Null secondary snapshot or node output pointer.");
    return OPENMC_E_INVALID_ARGUMENT;
  }
  const auto& nodes = static_cast<const SecondarySnapshot*>(snapshot)->nodes;
  if (index >= nodes.size()) {
    set_errmsg("Secondary snapshot node index is out of range.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  const auto& n = nodes[index];
  result->name = n.name.c_str();
  result->type = n.type;
  result->size = 0;
  result->data = nullptr;
  switch (n.type) {
  case OPENMC_SECONDARY_OBJECT:
  case OPENMC_SECONDARY_LIST:
    result->size = n.children.size();
    result->data = n.children.data();
    break;
  case OPENMC_SECONDARY_STRING:
    result->size = n.text.size();
    result->data = n.text.data();
    break;
  case OPENMC_SECONDARY_BOOL:
  case OPENMC_SECONDARY_INT:
  case OPENMC_SECONDARY_INT_ARRAY:
    result->size = n.integers.size();
    result->data = n.integers.data();
    break;
  case OPENMC_SECONDARY_DOUBLE:
  case OPENMC_SECONDARY_DOUBLE_ARRAY:
    result->size = n.reals.size();
    result->data = n.reals.data();
    break;
  }
  return 0;
}

extern "C" void openmc_secondary_snapshot_free(void* snapshot)
{
  delete static_cast<SecondarySnapshot*>(snapshot);
}
} // namespace openmc
