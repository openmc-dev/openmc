#include "openmc/geometry.h"

#include <array>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/lattice.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/string_utils.h"
#include "openmc/surface.h"


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================


namespace model {

int root_universe {-1};
int n_coord_levels;

std::vector<int64_t> overlap_check_count;

} // namespace model

//==============================================================================
// Non-member functions
//==============================================================================

bool check_cell_overlap(Particle* p, bool error)
{
  int n_coord = p->n_coord_;

  // Loop through each coordinate level
  for (int j = 0; j < n_coord; j++) {
    Universe& univ = *model::universes[p->coord_[j].universe];

    // Loop through each cell on this level
    for (auto index_cell : univ.cells_) {
      Cell& c = *model::cells[index_cell];
      if (c.contains(p->coord_[j].r, p->coord_[j].u, p->surface_)) {
        if (index_cell != p->coord_[j].cell) {
          if (error) {
            fatal_error(fmt::format(
              "Overlapping cells detected: {}, {} on universe {}",
              c.id_, model::cells[p->coord_[j].cell]->id_, univ.id_));
          }
          return true;
        }
        ++model::overlap_check_count[index_cell];
      }
    }
  }

  return false;
}

//==============================================================================

bool
find_cell_inner(Particle* p, const NeighborList* neighbor_list)
{
  // Find which cell of this universe the particle is in.  Use the neighbor list
  // to shorten the search if one was provided.
  bool found = false;
  int32_t i_cell;
  if (neighbor_list) {
    for (auto it = neighbor_list->cbegin(); it != neighbor_list->cend(); ++it) {
      i_cell = *it;

      // Make sure the search cell is in the same universe.
      int i_universe = p->coord_[p->n_coord_-1].universe;
      if (model::cells[i_cell]->universe_ != i_universe) continue;

      // Check if this cell contains the particle.
      Position r {p->r_local()};
      Direction u {p->u_local()};
      auto surf = p->surface_;
      if (model::cells[i_cell]->contains(r, u, surf)) {
        p->coord_[p->n_coord_-1].cell = i_cell;
        found = true;
        break;
      }
    }

  } else {
    int i_universe = p->coord_[p->n_coord_-1].universe;
    const auto& univ {*model::universes[i_universe]};
    const auto& cells {
      !univ.partitioner_
      ? model::universes[i_universe]->cells_
      : univ.partitioner_->get_cells(p->r_local(), p->u_local())
    };
    for (auto it = cells.cbegin(); it != cells.cend(); it++) {
      i_cell = *it;

      // Make sure the search cell is in the same universe.
      int i_universe = p->coord_[p->n_coord_-1].universe;
      if (model::cells[i_cell]->universe_ != i_universe) continue;

      // Check if this cell contains the particle.
      Position r {p->r_local()};
      Direction u {p->u_local()};
      auto surf = p->surface_;
      if (model::cells[i_cell]->contains(r, u, surf)) {
        p->coord_[p->n_coord_-1].cell = i_cell;
        found = true;
        break;
      }
    }
  }

  // Announce the cell that the particle is entering.
  if (found && (settings::verbosity >= 10 || p->trace_)) {
    auto msg = fmt::format("    Entering cell {}", model::cells[i_cell]->id_);
    write_message(msg, 1);
  }

  if (found) {
    Cell& c {*model::cells[i_cell]};
    if (c.type_ == Fill::MATERIAL) {
      //=======================================================================
      //! Found a material cell which means this is the lowest coord level.

      // Find the distribcell instance number.
      int offset = 0;
      if (c.distribcell_index_ >= 0) {
        for (int i = 0; i < p->n_coord_; i++) {
          const auto& c_i {*model::cells[p->coord_[i].cell]};
          if (c_i.type_ == Fill::UNIVERSE) {
            offset += c_i.offset_[c.distribcell_index_];
          } else if (c_i.type_ == Fill::LATTICE) {
            auto& lat {*model::lattices[p->coord_[i+1].lattice]};
            int i_xyz[3] {p->coord_[i+1].lattice_x,
                          p->coord_[i+1].lattice_y,
                          p->coord_[i+1].lattice_z};
            if (lat.are_valid_indices(i_xyz)) {
              offset += lat.offset(c.distribcell_index_, i_xyz);
            }
          }
        }
      }
      p->cell_instance_ = offset;

      // Set the material and temperature.
      p->material_last_ = p->material_;
      if (c.material_.size() > 1) {
        p->material_ = c.material_[p->cell_instance_];
      } else {
        p->material_ = c.material_[0];
      }
      p->sqrtkT_last_ = p->sqrtkT_;
      if (c.sqrtkT_.size() > 1) {
        p->sqrtkT_ = c.sqrtkT_[p->cell_instance_];
      } else {
        p->sqrtkT_ = c.sqrtkT_[0];
      }

      return true;

    } else if (c.type_ == Fill::UNIVERSE) {
      //========================================================================
      //! Found a lower universe, update this coord level then search the next.

      // Set the lower coordinate level universe.
      auto& coord {p->coord_[p->n_coord_]};
      coord.universe = c.fill_;

      // Set the position and direction.
      coord.r = p->r_local();
      coord.u = p->u_local();

      // Apply translation.
      coord.r -= c.translation_;

      // Apply rotation.
      if (!c.rotation_.empty()) {
        coord.rotate(c.rotation_);
      }

      // Update the coordinate level and recurse.
      ++p->n_coord_;
      return find_cell_inner(p, nullptr);

    } else if (c.type_ == Fill::LATTICE) {
      //========================================================================
      //! Found a lower lattice, update this coord level then search the next.

      Lattice& lat {*model::lattices[c.fill_]};

      // Set the position and direction.
      auto& coord {p->coord_[p->n_coord_]};
      coord.r = p->r_local();
      coord.u = p->u_local();

      // Apply translation.
      coord.r -= c.translation_;

      // Apply rotation.
      if (!c.rotation_.empty()) {
        coord.rotate(c.rotation_);
      }

      // Determine lattice indices.
      auto i_xyz = lat.get_indices(coord.r, coord.u);

      // Get local position in appropriate lattice cell
      coord.r = lat.get_local_position(coord.r, i_xyz);

      // Set lattice indices.
      coord.lattice = c.fill_;
      coord.lattice_x = i_xyz[0];
      coord.lattice_y = i_xyz[1];
      coord.lattice_z = i_xyz[2];

      // Set the lower coordinate level universe.
      if (lat.are_valid_indices(i_xyz)) {
        coord.universe = lat[i_xyz];
      } else {
        if (lat.outer_ != NO_OUTER_UNIVERSE) {
          coord.universe = lat.outer_;
        } else {
          warning(fmt::format("Particle {} is outside lattice {} but the "
            "lattice has no defined outer universe.", p->id_, lat.id_));
          return false;
        }
      }

      // Update the coordinate level and recurse.
      ++p->n_coord_;
      return find_cell_inner(p, nullptr);
    }
  }

  return found;
}

//==============================================================================

bool
find_cell(Particle* p, bool use_neighbor_lists)
{
  // Determine universe (if not yet set, use root universe).
  int i_universe = p->coord_[p->n_coord_-1].universe;
  if (i_universe == C_NONE) {
    p->coord_[0].universe = model::root_universe;
    p->n_coord_ = 1;
    i_universe = model::root_universe;
  }

  // Reset all the deeper coordinate levels.
  for (int i = p->n_coord_; i < p->coord_.size(); i++) {
    p->coord_[i].reset();
  }

  if (use_neighbor_lists) {
    // Get the cell this particle was in previously.
    auto coord_lvl = p->n_coord_ - 1;
    auto i_cell = p->coord_[coord_lvl].cell;
    Cell& c {*model::cells[i_cell]};

    // Search for the particle in that cell's neighbor list.  Return if we
    // found the particle.
    bool found = find_cell_inner(p, &c.neighbors_);
    if (found) return found;

    // The particle could not be found in the neighbor list.  Try searching all
    // cells in this universe, and update the neighbor list if we find a new
    // neighboring cell.
    found = find_cell_inner(p, nullptr);
    if (found) c.neighbors_.push_back(p->coord_[coord_lvl].cell);
    return found;

  } else {
    // Search all cells in this universe for the particle.
    return find_cell_inner(p, nullptr);
  }
}

//==============================================================================

void
cross_lattice(Particle* p, const BoundaryInfo& boundary)
{
  auto& coord {p->coord_[p->n_coord_ - 1]};
  auto& lat {*model::lattices[coord.lattice]};

  if (settings::verbosity >= 10 || p->trace_) {
    write_message(fmt::format(
      "    Crossing lattice {}. Current position ({},{},{}). r={}",
      lat.id_, coord.lattice_x, coord.lattice_y, coord.lattice_z, p->r()), 1);
  }

  // Set the lattice indices.
  coord.lattice_x += boundary.lattice_translation[0];
  coord.lattice_y += boundary.lattice_translation[1];
  coord.lattice_z += boundary.lattice_translation[2];
  std::array<int, 3> i_xyz {coord.lattice_x, coord.lattice_y, coord.lattice_z};

  // Set the new coordinate position.
  const auto& upper_coord {p->coord_[p->n_coord_ - 2]};
  const auto& cell {model::cells[upper_coord.cell]};
  Position r = upper_coord.r;
  r -= cell->translation_;
  if (!cell->rotation_.empty()) {
    r = r.rotate(cell->rotation_);
  }
  p->r_local() = lat.get_local_position(r, i_xyz);

  if (!lat.are_valid_indices(i_xyz)) {
    // The particle is outside the lattice.  Search for it from the base coords.
    p->n_coord_ = 1;
    bool found = find_cell(p, 0);
    if (!found && p->alive_) {
      p->mark_as_lost(fmt::format("Could not locate particle {} after "
        "crossing a lattice boundary", p->id_));
    }

  } else {
    // Find cell in next lattice element.
    p->coord_[p->n_coord_-1].universe = lat[i_xyz];
    bool found = find_cell(p, 0);

    if (!found) {
      // A particle crossing the corner of a lattice tile may not be found.  In
      // this case, search for it from the base coords.
      p->n_coord_ = 1;
      bool found = find_cell(p, 0);
      if (!found && p->alive_) {
        p->mark_as_lost(fmt::format("Could not locate particle {} after "
          "crossing a lattice boundary", p->id_));
      }
    }
  }
}

//==============================================================================

BoundaryInfo distance_to_boundary(Particle* p)
{
  BoundaryInfo info;
  double d_lat = INFINITY;
  double d_surf = INFINITY;
  int32_t level_surf_cross;
  std::array<int, 3> level_lat_trans {};

  // Loop over each coordinate level.
  for (int i = 0; i < p->n_coord_; i++) {
    const auto& coord {p->coord_[i]};
    Position r {coord.r};
    Direction u {coord.u};
    Cell& c {*model::cells[coord.cell]};

    // Find the oncoming surface in this cell and the distance to it.
    auto surface_distance = c.distance(r, u, p->surface_, p);
    d_surf = surface_distance.first;
    level_surf_cross = surface_distance.second;

    // Find the distance to the next lattice tile crossing.
    if (coord.lattice != C_NONE) {
      auto& lat {*model::lattices[coord.lattice]};
      std::array<int, 3> i_xyz {coord.lattice_x, coord.lattice_y, coord.lattice_z};
      //TODO: refactor so both lattice use the same position argument (which
      //also means the lat.type attribute can be removed)
      std::pair<double, std::array<int, 3>> lattice_distance;
      switch (lat.type_) {
        case LatticeType::rect:
          lattice_distance = lat.distance(r, u, i_xyz);
          break;
        case LatticeType::hex:
          auto& cell_above {model::cells[p->coord_[i-1].cell]};
          Position r_hex {p->coord_[i-1].r};
          r_hex -= cell_above->translation_;
          if (coord.rotated) {
            r_hex = r_hex.rotate(cell_above->rotation_);
          }
          r_hex.z = coord.r.z;
          lattice_distance = lat.distance(r_hex, u, i_xyz);
          break;
      }
      d_lat = lattice_distance.first;
      level_lat_trans = lattice_distance.second;

      if (d_lat < 0) {
        p->mark_as_lost(fmt::format(
          "Particle {} had a negative distance to a lattice boundary", p->id_));
      }
    }

    // If the boundary on this coordinate level is coincident with a boundary on
    // a higher level then we need to make sure that the higher level boundary
    // is selected.  This logic must consider floating point precision.
    double& d = info.distance;
    if (d_surf < d_lat - FP_COINCIDENT) {
      if (d == INFINITY || (d - d_surf)/d >= FP_REL_PRECISION) {
        d = d_surf;

        // If the cell is not simple, it is possible that both the negative and
        // positive half-space were given in the region specification. Thus, we
        // have to explicitly check which half-space the particle would be
        // traveling into if the surface is crossed
        if (c.simple_) {
          info.surface_index = level_surf_cross;
        } else {
          Position r_hit = r + d_surf * u;
          Surface& surf {*model::surfaces[std::abs(level_surf_cross)-1]};
          Direction norm = surf.normal(r_hit);
          if (u.dot(norm) > 0) {
            info.surface_index = std::abs(level_surf_cross);
          } else {
            info.surface_index = -std::abs(level_surf_cross);
          }
        }

        info.lattice_translation[0] = 0;
        info.lattice_translation[1] = 0;
        info.lattice_translation[2] = 0;
        info.coord_level = i + 1;
      }
    } else {
      if (d == INFINITY || (d - d_lat)/d >= FP_REL_PRECISION) {
        d = d_lat;
        info.surface_index = 0;
        info.lattice_translation = level_lat_trans;
        info.coord_level = i + 1;
      }
    }
  }
  return info;
}

//==============================================================================
// C API
//==============================================================================

extern "C" int
openmc_find_cell(const double* xyz, int32_t* index, int32_t* instance)
{
  Particle p;

  p.r() = Position{xyz};
  p.u() = {0.0, 0.0, 1.0};

  if (!find_cell(&p, false)) {
    set_errmsg(fmt::format("Could not find cell at position {}.", p.r()));
    return OPENMC_E_GEOMETRY;
  }

  *index = p.coord_[p.n_coord_-1].cell;
  *instance = p.cell_instance_;
  return 0;
}

extern "C" int openmc_global_bounding_box(double* llc, double* urc) {
  auto bbox = model::universes.at(model::root_universe)->bounding_box();

  // set lower left corner values
  llc[0] = bbox.xmin;
  llc[1] = bbox.ymin;
  llc[2] = bbox.zmin;

  // set upper right corner values
  urc[0] = bbox.xmax;
  urc[1] = bbox.ymax;
  urc[2] = bbox.zmax;

  return 0;
}

} // namespace openmc
