#include "openmc/particle_type.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>

#include "openmc/string_utils.h"

namespace openmc {
namespace {

constexpr const char* ATOMIC_SYMBOL[] = {"", "H", "He", "Li", "Be", "B", "C",
  "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
  "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
  "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
  "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr",
  "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
  "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
  "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
  "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg",
  "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

constexpr int MAX_Z =
  static_cast<int>(sizeof(ATOMIC_SYMBOL) / sizeof(ATOMIC_SYMBOL[0])) - 1;

bool is_integer_string(const std::string& s)
{
  if (s.empty())
    return false;
  size_t i = 0;
  if (s[0] == '-' || s[0] == '+') {
    if (s.size() == 1)
      return false;
    i = 1;
  }
  for (; i < s.size(); ++i) {
    if (!std::isdigit(static_cast<unsigned char>(s[i])))
      return false;
  }
  return true;
}

int atomic_number_from_symbol(std::string_view symbol)
{
  for (int z = 1; z <= MAX_Z; ++z) {
    if (symbol == ATOMIC_SYMBOL[z]) {
      return z;
    }
  }
  return 0;
}

bool parse_gnds_nuclide(std::string_view name, int& Z, int& A, int& m)
{
  if (name.empty())
    return false;

  size_t pos = 0;
  if (!std::isupper(static_cast<unsigned char>(name[pos])))
    return false;

  std::string symbol;
  symbol += name[pos++];
  if (pos < name.size() &&
      std::islower(static_cast<unsigned char>(name[pos]))) {
    symbol += name[pos++];
  }

  if (pos >= name.size() ||
      !std::isdigit(static_cast<unsigned char>(name[pos]))) {
    return false;
  }

  size_t a_start = pos;
  while (
    pos < name.size() && std::isdigit(static_cast<unsigned char>(name[pos]))) {
    ++pos;
  }
  A = std::stoi(std::string {name.substr(a_start, pos - a_start)});
  if (A <= 0 || A > 999)
    return false;

  m = 0;
  if (pos < name.size()) {
    if (name[pos] != '_' || pos + 2 >= name.size() || name[pos + 1] != 'm') {
      return false;
    }
    pos += 2;
    size_t m_start = pos;
    while (pos < name.size() &&
           std::isdigit(static_cast<unsigned char>(name[pos]))) {
      ++pos;
    }
    if (m_start == pos)
      return false;
    m = std::stoi(std::string {name.substr(m_start, pos - m_start)});
    if (m < 0 || m > 9)
      return false;
  }

  if (pos != name.size())
    return false;

  Z = atomic_number_from_symbol(symbol);
  return Z != 0;
}

// Helper to convert nuclear PDG number to nuclide name
std::string nuclide_name_from_pdg(int32_t pdg)
{
  int32_t code = pdg;
  int m = code % 10;
  int A = (code / 10) % 1000;
  int Z = (code / 10000) % 1000;

  if (Z <= 0 || Z > MAX_Z || A <= 0 || A > 999) {
    throw std::invalid_argument {
      "Invalid nuclear PDG number: " + std::to_string(pdg)};
  }

  std::string name = ATOMIC_SYMBOL[Z] + std::to_string(A);
  if (m > 0) {
    name += "_m" + std::to_string(m);
  }
  return name;
}

} // namespace

//==============================================================================
// ParticleType member function implementations
//==============================================================================

ParticleType::ParticleType(std::string_view str)
{
  std::string s {str};
  strtrim(s);
  if (s.empty()) {
    throw std::invalid_argument {"Particle string is empty."};
  }

  std::string lower = s;
  to_lower(lower);

  // Check for pdg: prefix
  if (starts_with(lower, "pdg:")) {
    std::string value_str = lower.substr(4);
    if (!is_integer_string(value_str)) {
      throw std::invalid_argument {"Invalid PDG number: " + value_str};
    }
    pdg_number_ = std::stoi(value_str);
    return;
  }

  // Check for known particle names
  if (lower == "neutron" || lower == "n") {
    pdg_number_ = PDG_NEUTRON;
    return;
  }
  if (lower == "photon" || lower == "gamma") {
    pdg_number_ = PDG_PHOTON;
    return;
  }
  if (lower == "electron") {
    pdg_number_ = PDG_ELECTRON;
    return;
  }
  if (lower == "positron") {
    pdg_number_ = PDG_POSITRON;
    return;
  }
  if (lower == "proton" || lower == "p" || lower == "h1") {
    pdg_number_ = PDG_PROTON;
    return;
  }
  if (lower == "deuteron" || lower == "d" || lower == "h2") {
    pdg_number_ = PDG_DEUTERON;
    return;
  }
  if (lower == "triton" || lower == "t" || lower == "h3") {
    pdg_number_ = PDG_TRITON;
    return;
  }
  if (lower == "alpha" || lower == "he4") {
    pdg_number_ = PDG_ALPHA;
    return;
  }

  // Check for integer string
  if (is_integer_string(s)) {
    pdg_number_ = std::stoi(s);
    return;
  }

  // Try to parse as GNDS nuclide name
  int Z = 0;
  int A = 0;
  int m = 0;
  if (!parse_gnds_nuclide(s, Z, A, m)) {
    throw std::invalid_argument {"Invalid nuclide name: " + s};
  }
  pdg_number_ = 1000000000 + Z * 10000 + A * 10 + m;
}

std::string ParticleType::str() const
{
  if (pdg_number_ == PDG_NEUTRON)
    return "neutron";
  if (pdg_number_ == PDG_PHOTON)
    return "photon";
  if (pdg_number_ == PDG_ELECTRON)
    return "electron";
  if (pdg_number_ == PDG_POSITRON)
    return "positron";
  if (pdg_number_ == PDG_PROTON)
    return "proton";

  if (is_nucleus()) {
    return nuclide_name_from_pdg(pdg_number_);
  }

  return "pdg:" + std::to_string(pdg_number_);
}

//==============================================================================
// Free function implementations
//==============================================================================

ParticleType legacy_particle_index_to_type(int index)
{
  switch (index) {
  case 0:
    return ParticleType {PDG_NEUTRON};
  case 1:
    return ParticleType {PDG_PHOTON};
  case 2:
    return ParticleType {PDG_ELECTRON};
  case 3:
    return ParticleType {PDG_POSITRON};
  default:
    throw std::invalid_argument {
      "Invalid legacy particle index: " + std::to_string(index)};
  }
}

} // namespace openmc
