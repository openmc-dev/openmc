#include <algorithm> // for min, max, sort, fill
#include <cassert>
#include <cmath>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_set>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xoperation.hpp"
#include "xtensor/xview.hpp"

#include "openmc/capi.h"
#include "openmc/container_util.h"
#include "openmc/cross_sections.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/math_functions.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/stochastic_media.h"
#include "openmc/string_utils.h"
#include "openmc/thermal.h"
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
    if (mats.size() > 0) {
      matrix_mat_.reserve(mats.size());
      for (std::string mat : mats) {
        if (mat.compare("void") == 0) {
          matrix_mat_.push_back(MATERIAL_VOID);
        } else {
          matrix_mat_.push_back(std::stoi(mat));
        }
      }
    }
  } else {
    fatal_error(fmt::format(
      "An empty material element was specified for stochastic media {}", id_));
  }

  if (check_for_node(node, "radius")) {
    radius_ = std::stoi(get_node_value(node, "radius"));
  } else {
    fatal_error(fmt::format(
      "An empty particle radius was specified for stochastic media {}", id_));
  }

  if (check_for_node(node, "pack_fraction")) {
    pf_ = std::stoi(get_node_value(node, "pack_fraction"));
  } else {
    fatal_error(fmt::format(
      "An empty pack fraction was specified for stochastic media {}", id_));
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