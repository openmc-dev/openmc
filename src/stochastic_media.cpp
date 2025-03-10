#include <algorithm> // for min, max, sort, fill
#include <cassert>
#include <cmath>
#include <iterator>

#include "openmc/cell.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/material.h"
#include "openmc/math_functions.h"
#include "openmc/simulation.h"
#include "openmc/stochastic_media.h"
#include "openmc/string_utils.h"
#include "openmc/xml_interface.h"

namespace openmc {
//==============================================================================
// Global variables
//==============================================================================
namespace model {
std::unordered_map<int32_t, int32_t> stochastic_media_map;
vector<unique_ptr<Stochastic_Media>> stochastic_media;

} // namespace model
//==============================================================================
// Stochastic_Media implementation
//==============================================================================

Stochastic_Media::~Stochastic_Media()
{
  model::stochastic_media_map.erase(id_);
}

void Stochastic_Media::adjust_indices()
{
  // Adjust the indices for the materials array.
  for (auto& mat : particle_mat_) {

    auto search = model::material_map.find(mat);
    if (search != model::universe_map.end()) {
      mat = search->second;
    } else {
      fatal_error(fmt::format("Invalid material number {} specified on the "
                              "particle for stochastic media {}",
        mat, id_));
    }
  }

  auto search = model::material_map.find(matrix_mat_);
  if (search != model::universe_map.end()) {
    matrix_mat_ = search->second;
  } else {
    fatal_error(fmt::format("Invalid material number {} specified on the "
                            "matrix for stochastic media {}",
      matrix_mat_, id_));
  }
}

CLS_Media::CLS_Media(pugi::xml_node node)
{
  index_ = model::stochastic_media.size(); // Avoids warning about narrowing
  if (check_for_node(node, "id")) {
    this->set_id(std::stoi(get_node_value(node, "id")));
  } else {
    fatal_error("Must specify id of material in materials XML file.");
  }

  if (check_for_node(node, "name")) {
    name_ = get_node_value(node, "name");
  }

  bool particle_mat_present = check_for_node(node, "particle_material");
  bool matrix_mat_present = check_for_node(node, "matrix_material");

  // Read the material element.  There can be zero materials (filled with a
  // universe), more than one material (distribmats), and some materials may
  // be "void".
  if (particle_mat_present) {
    vector<std::string> mats {
      get_node_array<std::string>(node, "particle_material", true)};
    if (mats.size() > 0) {
      particle_mat_.reserve(mats.size());
      for (std::string mat : mats) {
        if (mat.compare("void") == 0) {
          particle_mat_.push_back(MATERIAL_VOID);
        } else {
          particle_mat_.push_back(std::stoi(mat));
        }
      }
    }
  } else {
    fatal_error(fmt::format(
      "An empty material element was specified for stochastic media {}", id_));
  }

  if (matrix_mat_present) {
    vector<std::string> mats {
      get_node_array<std::string>(node, "matrix_material", true)};
    matrix_mat_ = std::stoi(mats[0]);
    if (mats.size() > 1) {
      fatal_error(fmt::format(
        "Only one matrix material can be specified for stochastic media {}",
        id_));
    }
  } else {
    fatal_error(fmt::format(
      "An empty matrix material was specified for stochastic media {}", id_));
  }

  if (check_for_node(node, "radius")) {
    radius_ = node.attribute("radius").as_double();
  } else {
    fatal_error(fmt::format(
      "An empty particle radius was specified for stochastic media {}", id_));
  }

  if (check_for_node(node, "pack_fraction")) {
    pf_ = node.attribute("pack_fraction").as_double();
  } else {
    fatal_error(fmt::format(
      "An empty pack fraction was specified for stochastic media {}", id_));
  }
}

double CLS_Media::distance_to_stochamedia(Particle& p)
{ // Sample the distance to the stochastic media
  double distance = INFINITY;
  if (p.status() == ParticleStatus::IN_STOCHASTIC_MEDIA) {
    // Designed for a randomized medium only for the time
    // being, to be upgraded subsequently
    double cos_value = sqrt(prn(p.current_seed()));
    distance = 2 * radius_ * cos_value;
    p.status() = ParticleStatus::IN_STOCHASTIC_MEDIA;
  } else if (p.status() == ParticleStatus::IN_MATRIX) {
    double matrix_mean_chord = 4 / 3 * radius_ * (1 - pf_) / pf_;
    distance = -matrix_mean_chord * std::log(prn(p.current_seed()));
  }
  return distance;
}
void CLS_Media::sample_material(Particle& p)
{
  p.material_last() = p.material();
  p.sqrtkT_last() = p.sqrtkT();

  // Sample the material based on the packing fraction
  auto i_cell = p.lowest_coord().cell;
  Cell& c {*model::cells[i_cell]};
  if (p.status() == ParticleStatus::OUTSIDE) {
    double rand = openmc::prn(p.current_seed());
    if (rand < pf_) {
      p.status() = ParticleStatus::IN_STOCHASTIC_MEDIA;

      p.material() = this->particle_mat(0);

      p.sqrtkT() = c.sqrtkT(p.cell_instance());
    } else {
      p.status() = ParticleStatus::IN_MATRIX;
    }
  }
}

void Stochastic_Media::set_id(int32_t id)
{
  assert(id >= 0 || id == C_NONE);

  // Clear entry in stochastic media map if an ID was already assigned before
  if (id_ != C_NONE) {
    model::stochastic_media_map.erase(id_);
    id_ = C_NONE;
  }

  // Make sure no other material has same ID
  if (model::stochastic_media_map.find(id) !=
      model::stochastic_media_map.end()) {
    throw std::runtime_error {
      "Two stochastic media have the same ID: " + std::to_string(id)};
  }

  // If no ID specified, auto-assign next ID in sequence
  if (id == C_NONE) {
    id = 0;
    for (const auto& m : model::stochastic_media) {
      id = std::max(id, m->id_);
    }
    ++id;
  }

  // Update ID and entry in material map
  id_ = id;
  model::stochastic_media_map[id] = index_;
}
//==============================================================================
// Non-method functions
//==============================================================================
void read_stochastic_media(pugi::xml_node root)
{
  // read stochastic media
  if (check_for_node(root, "CLS_media")) {
    for (pugi::xml_node media_node : root.children("CLS_media")) {
      model::stochastic_media.push_back(make_unique<CLS_Media>(media_node));
    }
    model::stochastic_media.shrink_to_fit();
  }
}

void free_memory_stochastic_media()
{
  model::stochastic_media.clear();
  model::stochastic_media_map.clear();
}
} // namespace openmc