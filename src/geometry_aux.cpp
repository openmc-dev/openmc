#include "openmc/geometry_aux.h"

#include <algorithm> // for std::max
#include <sstream>
#include <unordered_set>

#include <fmt/core.h>
#include <pugixml.hpp>

#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/container_util.h"
#include "openmc/dagmc.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/lattice.h"
#include "openmc/material.h"
#include "openmc/settings.h"
#include "openmc/surface.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_cell_instance.h"
#include "openmc/tallies/filter_distribcell.h"

namespace openmc {

namespace model {
std::unordered_map<int32_t, int32_t> universe_level_counts;
} // namespace model

void read_geometry_xml()
{
  // Display output message
  write_message("Reading geometry XML file...", 5);

  // Check if geometry.xml exists
  std::string filename = settings::path_input + "geometry.xml";
  if (!file_exists(filename)) {
    fatal_error("Geometry XML file '" + filename + "' does not exist!");
  }

  // Parse settings.xml file
  pugi::xml_document doc;
  auto result = doc.load_file(filename.c_str());
  if (!result) {
    fatal_error("Error processing geometry.xml file.");
  }

  // Get root element
  pugi::xml_node root = doc.document_element();

  read_geometry_xml(root);
}

void read_geometry_xml(pugi::xml_node root)
{
  // Read surfaces, cells, lattice
  std::set<std::pair<int, int>> periodic_pairs;
  std::unordered_map<int, double> albedo_map;
  std::unordered_map<int, int> periodic_sense_map;

  read_surfaces(root, periodic_pairs, albedo_map, periodic_sense_map);
  read_cells(root);
  prepare_boundary_conditions(periodic_pairs, albedo_map, periodic_sense_map);
  read_lattices(root);

  // Check to make sure a boundary condition was applied to at least one
  // surface
  bool boundary_exists = false;
  for (const auto& surf : model::surfaces) {
    if (surf->bc_) {
      boundary_exists = true;
      break;
    }
  }

  if (settings::run_mode != RunMode::PLOTTING &&
      settings::run_mode != RunMode::VOLUME && !boundary_exists) {
    fatal_error("No boundary conditions were applied to any surfaces!");
  }

  // Allocate universes, universe cell arrays, and assign base universe
  model::root_universe = find_root_universe();

  // if the root universe is DAGMC geometry, make sure the model is well-formed
  check_dagmc_root_univ();
}

//==============================================================================

void adjust_indices()
{
  // Adjust material/fill idices.
  for (auto& c : model::cells) {
    if (c->fill_ != C_NONE) {
      int32_t id = c->fill_;
      auto search_univ = model::universe_map.find(id);
      auto search_lat = model::lattice_map.find(id);
      if (search_univ != model::universe_map.end()) {
        c->type_ = Fill::UNIVERSE;
        c->fill_ = search_univ->second;
      } else if (search_lat != model::lattice_map.end()) {
        c->type_ = Fill::LATTICE;
        c->fill_ = search_lat->second;
      } else {
        fatal_error(fmt::format("Specified fill {} on cell {} is neither a "
                                "universe nor a lattice.",
          id, c->id_));
      }
    } else {
      c->type_ = Fill::MATERIAL;
      for (auto& mat_id : c->material_) {
        if (mat_id != MATERIAL_VOID) {
          auto search = model::material_map.find(mat_id);
          if (search == model::material_map.end()) {
            fatal_error(
              fmt::format("Could not find material {} specified on cell {}",
                mat_id, c->id_));
          }
          // Change from ID to index
          mat_id = search->second;
        }
      }
    }
  }

  // Change cell.universe values from IDs to indices.
  for (auto& c : model::cells) {
    auto search = model::universe_map.find(c->universe_);
    if (search != model::universe_map.end()) {
      c->universe_ = search->second;
    } else {
      fatal_error(fmt::format("Could not find universe {} specified on cell {}",
        c->universe_, c->id_));
    }
  }

  // Change all lattice universe values from IDs to indices.
  for (auto& l : model::lattices) {
    l->adjust_indices();
  }
}

//==============================================================================
//! Partition some universes with many z-planes for faster find_cell searches.

void partition_universes()
{
  // Iterate over universes with more than 10 cells.  (Fewer than 10 is likely
  // not worth partitioning.)
  for (const auto& univ : model::universes) {
    if (univ->cells_.size() > 10) {
      // Collect the set of surfaces in this universe.
      std::unordered_set<int32_t> surf_inds;
      for (auto i_cell : univ->cells_) {
        for (auto token : model::cells[i_cell]->surfaces()) {
          surf_inds.insert(std::abs(token) - 1);
        }
      }

      // Partition the universe if there are more than 5 z-planes.  (Fewer than
      // 5 is likely not worth it.)
      int n_zplanes = 0;
      for (auto i_surf : surf_inds) {
        if (dynamic_cast<const SurfaceZPlane*>(model::surfaces[i_surf].get())) {
          ++n_zplanes;
          if (n_zplanes > 5) {
            univ->partitioner_ = make_unique<UniversePartitioner>(*univ);
            break;
          }
        }
      }
    }
  }
}

//==============================================================================

void assign_temperatures()
{
  for (auto& c : model::cells) {
    // Ignore non-material cells and cells with defined temperature.
    if (c->material_.size() == 0)
      continue;
    if (c->sqrtkT_.size() > 0)
      continue;

    c->sqrtkT_.reserve(c->material_.size());
    for (auto i_mat : c->material_) {
      if (i_mat == MATERIAL_VOID) {
        // Set void region to 0K.
        c->sqrtkT_.push_back(0);
      } else {
        const auto& mat {model::materials[i_mat]};
        c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * mat->temperature()));
      }
    }
  }
}

//==============================================================================

void finalize_cell_densities()
{
  for (auto& c : model::cells) {
    // Convert to density multipliers.
    if (!c->density_mult_.empty()) {
      for (int32_t instance = 0; instance < c->density_mult_.size();
           ++instance) {
        c->density_mult_[instance] /=
          model::materials[c->material(instance)]->density_gpcc();
      }
    } else {
      c->density_mult_ = {1.0};
    }
  }
}

//==============================================================================

void get_temperatures(
  vector<vector<double>>& nuc_temps, vector<vector<double>>& thermal_temps)
{
  for (const auto& cell : model::cells) {
    // Skip non-material cells.
    if (cell->fill_ != C_NONE)
      continue;

    for (int j = 0; j < cell->material_.size(); ++j) {
      // Skip void materials
      int i_material = cell->material_[j];
      if (i_material == MATERIAL_VOID)
        continue;

      // Get temperature(s) of cell (rounding to nearest integer)
      vector<double> cell_temps;
      if (cell->sqrtkT_.size() == 1) {
        double sqrtkT = cell->sqrtkT_[0];
        cell_temps.push_back(sqrtkT * sqrtkT / K_BOLTZMANN);
      } else if (cell->sqrtkT_.size() == cell->material_.size()) {
        double sqrtkT = cell->sqrtkT_[j];
        cell_temps.push_back(sqrtkT * sqrtkT / K_BOLTZMANN);
      } else {
        for (double sqrtkT : cell->sqrtkT_)
          cell_temps.push_back(sqrtkT * sqrtkT / K_BOLTZMANN);
      }

      const auto& mat {model::materials[i_material]};
      for (const auto& i_nuc : mat->nuclide_) {
        for (double temperature : cell_temps) {
          // Add temperature if it hasn't already been added
          if (!contains(nuc_temps[i_nuc], temperature))
            nuc_temps[i_nuc].push_back(temperature);
        }
      }

      for (const auto& table : mat->thermal_tables_) {
        // Get index in data::thermal_scatt array
        int i_sab = table.index_table;

        for (double temperature : cell_temps) {
          // Add temperature if it hasn't already been added
          if (!contains(thermal_temps[i_sab], temperature))
            thermal_temps[i_sab].push_back(temperature);
        }
      }
    }
  }
}

//==============================================================================

void finalize_geometry()
{
  // Perform some final operations to set up the geometry
  adjust_indices();
  count_universe_instances();
  partition_universes();

  // Assign temperatures to cells that don't have temperatures already assigned
  assign_temperatures();

  // Determine number of nested coordinate levels in the geometry
  model::n_coord_levels = maximum_levels(model::root_universe);
}

//==============================================================================

int32_t find_root_universe()
{
  // Find all the universes listed as a cell fill.
  std::unordered_set<int32_t> fill_univ_ids;
  for (const auto& c : model::cells) {
    fill_univ_ids.insert(c->fill_);
  }

  // Find all the universes contained in a lattice.
  for (const auto& lat : model::lattices) {
    for (auto it = lat->begin(); it != lat->end(); ++it) {
      fill_univ_ids.insert(*it);
    }
    if (lat->outer_ != NO_OUTER_UNIVERSE) {
      fill_univ_ids.insert(lat->outer_);
    }
  }

  // Figure out which universe is not in the set.  This is the root universe.
  bool root_found {false};
  int32_t root_univ;
  for (int32_t i = 0; i < model::universes.size(); i++) {
    auto search = fill_univ_ids.find(model::universes[i]->id_);
    if (search == fill_univ_ids.end()) {
      if (root_found) {
        fatal_error("Two or more universes are not used as fill universes, so "
                    "it is not possible to distinguish which one is the root "
                    "universe.");
      } else {
        root_found = true;
        root_univ = i;
      }
    }
  }
  if (!root_found)
    fatal_error("Could not find a root universe.  Make sure "
                "there are no circular dependencies in the geometry.");

  return root_univ;
}

//==============================================================================

DistribcellOffsets::DistribcellOffsets(
  const vector<int32_t>& maps, const int32_t* values)
  : values_(values), maps_(&maps)
{
  if (!maps.empty()) {
    first_map_ = maps[0];
    n_first_ = 1;
    while (n_first_ < maps.size() && maps[n_first_] == first_map_ + n_first_) {
      ++n_first_;
    }
  }
}

int32_t DistribcellOffsets::find(int32_t map) const
{
  if (!maps_)
    return 0;
  auto it = std::lower_bound(maps_->begin(), maps_->end(), map);
  if (it == maps_->end() || *it != map)
    return 0;
  return values_[it - maps_->begin()];
}

bool DistribcellOffsets::contains(int32_t map) const
{
  return maps_ && std::binary_search(maps_->begin(), maps_->end(), map);
}

//==============================================================================

namespace {

//! Builds the distributed cell offsets, visiting each universe and lattice
//! once, after the universes they contain.
class DistribcellBuilder {
public:
  explicit DistribcellBuilder(const vector<bool>& is_target)
    : is_target_(is_target), univ_map_(model::universes.size(), C_NONE),
      univ_visited_(model::universes.size(), false),
      lat_visited_(model::lattices.size(), false),
      univ_counts_(model::universes.size()), lat_counts_(model::lattices.size())
  {}

  //! Visit the geometry from the root universe. Maps are numbered in
  //! post-order, so that the maps in a universe are mostly consecutive.
  void build() { visit_universe(model::root_universe); }

  int32_t map(int32_t univ) const { return univ_map_[univ]; }

private:
  //! Find the maps in a universe, store the offsets of its fill cells, and
  //! count the instances of each map.
  void visit_universe(int32_t univ_indx)
  {
    if (univ_visited_[univ_indx])
      return;
    univ_visited_[univ_indx] = true;
    Universe& univ = *model::universes[univ_indx];

    auto& maps = univ.distribcell_maps_;
    maps.clear();
    for (int32_t cell_indx : univ.cells_) {
      const Cell& c = *model::cells[cell_indx];
      if (c.type_ == Fill::UNIVERSE) {
        visit_universe(c.fill_);
      } else if (c.type_ == Fill::LATTICE) {
        visit_lattice(c.fill_);
      } else {
        continue;
      }
      const auto& fill_maps = this->fill_maps(c);
      maps.insert(maps.end(), fill_maps.begin(), fill_maps.end());
    }
    sort_unique(maps);

    // The universe's own map comes after the maps it contains
    if (is_target_[univ_indx]) {
      univ_map_[univ_indx] = n_maps_++;
      maps.push_back(univ_map_[univ_indx]);
    }

    // The instances counted before a cell are its offsets. The offsets of all
    // cells are stored together, so the storage is allocated once.
    size_t n_values = 0;
    for (int32_t cell_indx : univ.cells_) {
      const Cell& c = *model::cells[cell_indx];
      if (c.type_ != Fill::MATERIAL)
        n_values += fill_maps(c).size();
    }
    univ.offset_values_.clear();
    univ.offset_values_.reserve(n_values);
    auto& counts = univ_counts_[univ_indx];
    counts.assign(maps.size(), 0);
    for (int32_t cell_indx : univ.cells_) {
      Cell& c = *model::cells[cell_indx];
      if (c.type_ == Fill::MATERIAL)
        continue;
      const auto& fill_maps = this->fill_maps(c);
      c.offset_ = store(univ.offset_values_, maps, counts, fill_maps);
      add_counts(maps, counts, fill_maps,
        c.type_ == Fill::UNIVERSE ? univ_counts_[c.fill_]
                                  : lat_counts_[c.fill_]);
    }
    if (is_target_[univ_indx])
      counts.back() = 1;
  }

  //! Find the maps in a lattice, store the offsets of its tiles, and count the
  //! instances of each map. The outer universe is not counted, but its maps are
  //! included since a lattice cell's offsets are read for particles in the
  //! outer universe.
  void visit_lattice(int32_t lat_indx)
  {
    if (lat_visited_[lat_indx])
      return;
    lat_visited_[lat_indx] = true;
    Lattice& lat = *model::lattices[lat_indx];

    auto& maps = lat.distribcell_maps_;
    maps.clear();
    std::unordered_set<int32_t> univs;
    for (LatticeIter it = lat.begin(); it != lat.end(); ++it) {
      if (univs.insert(*it).second) {
        visit_universe(*it);
        const auto& univ_maps = model::universes[*it]->distribcell_maps_;
        maps.insert(maps.end(), univ_maps.begin(), univ_maps.end());
      }
    }
    if (lat.outer_ != NO_OUTER_UNIVERSE) {
      visit_universe(lat.outer_);
      const auto& outer_maps = model::universes[lat.outer_]->distribcell_maps_;
      maps.insert(maps.end(), outer_maps.begin(), outer_maps.end());
    }
    sort_unique(maps);

    // The instances counted before a tile are its offsets. The offsets of all
    // tiles are stored together, so the storage is allocated once.
    size_t n_values = 0;
    for (LatticeIter it = lat.begin(); it != lat.end(); ++it) {
      n_values += model::universes[*it]->distribcell_maps_.size();
    }
    lat.offset_values_.clear();
    lat.offset_values_.reserve(n_values);
    auto& counts = lat_counts_[lat_indx];
    counts.assign(maps.size(), 0);
    lat.offsets_.assign(lat.universes_.size(), {});
    for (LatticeIter it = lat.begin(); it != lat.end(); ++it) {
      const auto& univ_maps = model::universes[*it]->distribcell_maps_;
      lat.offsets_[it.indx_] =
        store(lat.offset_values_, maps, counts, univ_maps);
      add_counts(maps, counts, univ_maps, univ_counts_[*it]);
    }
  }

  //! Sorted maps in the fill of a fill cell
  static const vector<int32_t>& fill_maps(const Cell& c)
  {
    return c.type_ == Fill::UNIVERSE
             ? model::universes[c.fill_]->distribcell_maps_
             : model::lattices[c.fill_]->distribcell_maps_;
  }

  //! Append the counts of a subset of maps to storage with enough capacity
  //! \return Offsets of the subset
  static DistribcellOffsets store(vector<int32_t>& storage,
    const vector<int32_t>& maps, const vector<int32_t>& counts,
    const vector<int32_t>& subset)
  {
    const int32_t* values = storage.data() + storage.size();
    for (int32_t m : subset) {
      storage.push_back(counts[position(maps, m)]);
    }
    return DistribcellOffsets(subset, values);
  }

  //! Add the counts of a subset of maps to counts
  static void add_counts(const vector<int32_t>& maps, vector<int32_t>& counts,
    const vector<int32_t>& subset, const vector<int32_t>& subset_counts)
  {
    for (int32_t i = 0; i < subset.size(); ++i) {
      counts[position(maps, subset[i])] += subset_counts[i];
    }
  }

  //! Position of a map in a sorted list of maps that contains it
  static int32_t position(const vector<int32_t>& maps, int32_t map)
  {
    return std::lower_bound(maps.begin(), maps.end(), map) - maps.begin();
  }

  static void sort_unique(vector<int32_t>& v)
  {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  }

  const vector<bool>& is_target_;
  int32_t n_maps_ {0};
  vector<int32_t> univ_map_;
  vector<bool> univ_visited_;
  vector<bool> lat_visited_;
  // Number of instances of each map in a universe or lattice
  vector<vector<int32_t>> univ_counts_;
  vector<vector<int32_t>> lat_counts_;
};

} // namespace

void prepare_distribcell(const std::vector<int32_t>* user_distribcells)
{
  write_message("Preparing distributed cell instances...", 5);

  std::unordered_set<int32_t> distribcells;

  // start with any cells manually specified via the C++ API
  if (user_distribcells) {
    distribcells.insert(user_distribcells->begin(), user_distribcells->end());
  }

  // Find all cells listed in a DistribcellFilter or CellInstanceFilter
  for (auto& filt : model::tally_filters) {
    auto* distrib_filt = dynamic_cast<DistribcellFilter*>(filt.get());
    auto* cell_inst_filt = dynamic_cast<CellInstanceFilter*>(filt.get());
    if (distrib_filt) {
      distribcells.insert(distrib_filt->cell());
    }
    if (cell_inst_filt) {
      const auto& filter_cells = cell_inst_filt->cells();
      distribcells.insert(filter_cells.begin(), filter_cells.end());
    }
  }

  // By default, add material cells to the list of distributed cells
  if (settings::material_cell_offsets) {
    for (int64_t i = 0; i < model::cells.size(); ++i) {
      if (model::cells[i]->type_ == Fill::MATERIAL)
        distribcells.insert(i);
    }
  }

  // Make sure that the number of materials/temperatures matches the number of
  // cell instances.
  for (int i = 0; i < model::cells.size(); i++) {
    Cell& c {*model::cells[i]};

    if (c.material_.size() > 1) {
      if (c.material_.size() != c.n_instances()) {
        fatal_error(fmt::format(
          "Cell {} was specified with {} materials but has {} distributed "
          "instances. The number of materials must equal one or the number "
          "of instances.",
          c.id_, c.material_.size(), c.n_instances()));
      }
    }

    if (c.sqrtkT_.size() > 1) {
      if (c.sqrtkT_.size() != c.n_instances()) {
        fatal_error(fmt::format(
          "Cell {} was specified with {} temperatures but has {} distributed "
          "instances. The number of temperatures must equal one or the number "
          "of instances.",
          c.id_, c.sqrtkT_.size(), c.n_instances()));
      }
    }

    if (c.density_mult_.size() > 1) {
      if (c.density_mult_.size() != c.n_instances()) {
        fatal_error(fmt::format("Cell {} was specified with {} density "
                                "multipliers but has {} distributed "
                                "instances. The number of density multipliers "
                                "must equal one or the number "
                                "of instances.",
          c.id_, c.density_mult_.size(), c.n_instances()));
      }
    }
  }

  // Each universe containing distributed cells gets a map, which numbers the
  // instances of that universe.
  vector<bool> is_target(model::universes.size(), false);
  for (auto idx : distribcells) {
    is_target[model::cells[idx]->universe_] = true;
  }

  DistribcellBuilder builder(is_target);
  builder.build();

  for (auto idx : distribcells) {
    Cell& c = *model::cells[idx];
    c.distribcell_index_ = builder.map(c.universe_);
  }
}

//==============================================================================

void count_universe_instances()
{
  // Call a function with each universe filling a cell or lattice element of a
  // universe, once per cell or lattice element
  auto for_each_fill = [](int32_t i_univ, auto&& f) {
    for (int32_t i_cell : model::universes[i_univ]->cells_) {
      Cell& c = *model::cells[i_cell];
      if (c.type_ == Fill::UNIVERSE) {
        f(c.fill_);
      } else if (c.type_ == Fill::LATTICE) {
        Lattice& lat = *model::lattices[c.fill_];
        for (auto it = lat.begin(); it != lat.end(); ++it)
          f(*it);
      }
    }
  };

  // Order the universes reachable from the root so that every universe comes
  // after all universes that contain it
  vector<int32_t> order;
  vector<bool> visited(model::universes.size(), false);
  auto visit = [&](int32_t i_univ, auto& self) -> void {
    visited[i_univ] = true;
    for_each_fill(i_univ, [&](int32_t next) {
      if (!visited[next])
        self(next, self);
    });
    order.push_back(i_univ);
  };
  visit(model::root_universe, visit);

  // The number of instances of a universe is the sum over the cells and
  // lattice elements it fills of the number of instances of the universe
  // containing them. Universes not reachable from the root have none.
  for (auto& univ : model::universes)
    univ->n_instances_ = 0;
  model::universes[model::root_universe]->n_instances_ = 1;
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    int n = model::universes[*it]->n_instances_;
    for_each_fill(
      *it, [&](int32_t next) { model::universes[next]->n_instances_ += n; });
  }
}

//==============================================================================

std::string distribcell_path_inner(int32_t target_cell, int32_t map,
  int32_t target_offset, const Universe& search_univ, int32_t offset)
{
  std::stringstream path;

  path << "u" << search_univ.id_ << "->";

  // Check to see if this universe directly contains the target cell.  If so,
  // write to the path and return.
  for (int32_t cell_indx : search_univ.cells_) {
    if ((cell_indx == target_cell) && (offset == target_offset)) {
      Cell& c = *model::cells[cell_indx];
      path << "c" << c.id_;
      return path.str();
    }
  }

  // The target must be further down the geometry tree and contained in a fill
  // cell or lattice cell in this universe.  Find which cell contains the
  // target.
  vector<std::int32_t>::const_reverse_iterator cell_it {
    search_univ.cells_.crbegin()};
  for (; cell_it != search_univ.cells_.crend(); ++cell_it) {
    Cell& c = *model::cells[*cell_it];

    // Material cells don't contain other cells and other cells may not contain
    // the target's universe, so ignore them.
    if (c.type_ != Fill::MATERIAL && c.offset_.contains(map)) {
      int32_t temp_offset = offset + c.offset_[map];

      // The desired cell is the first cell that gives an offset smaller or
      // equal to the target offset.
      if (temp_offset <= target_offset)
        break;
    }
  }

  // if we get through the loop without finding an appropriate entry, throw
  // an error
  if (cell_it == search_univ.cells_.crend()) {
    fatal_error(
      fmt::format("Failed to generate a text label for distribcell with ID {}."
                  "The current label is: '{}'",
        model::cells[target_cell]->id_, path.str()));
  }

  // Add the cell to the path string.
  Cell& c = *model::cells[*cell_it];
  path << "c" << c.id_ << "->";

  if (c.type_ == Fill::UNIVERSE) {
    // Recurse into the fill cell.
    offset += c.offset_[map];
    path << distribcell_path_inner(
      target_cell, map, target_offset, *model::universes[c.fill_], offset);
    return path.str();
  } else {
    // Recurse into the lattice cell.
    Lattice& lat = *model::lattices[c.fill_];
    path << "l" << lat.id_;
    for (ReverseLatticeIter it = lat.rbegin(); it != lat.rend(); ++it) {
      if (!lat.offsets_[it.indx_].contains(map))
        continue;
      int32_t temp_offset = offset + lat.offset(map, it.indx_) + c.offset_[map];
      if (temp_offset <= target_offset) {
        offset = temp_offset;
        path << "(" << lat.index_to_string(it.indx_) << ")->";
        path << distribcell_path_inner(
          target_cell, map, target_offset, *model::universes[*it], offset);
        return path.str();
      }
    }
    throw std::runtime_error {"Error determining distribcell path."};
  }
}

std::string distribcell_path(
  int32_t target_cell, int32_t map, int32_t target_offset)
{
  auto& root_univ = *model::universes[model::root_universe];
  return distribcell_path_inner(target_cell, map, target_offset, root_univ, 0);
}

//==============================================================================

int maximum_levels(int32_t univ)
{

  const auto level_count = model::universe_level_counts.find(univ);
  if (level_count != model::universe_level_counts.end()) {
    return level_count->second;
  }

  int levels_below {0};

  for (int32_t cell_indx : model::universes[univ]->cells_) {
    Cell& c = *model::cells[cell_indx];
    if (c.type_ == Fill::UNIVERSE) {
      int32_t next_univ = c.fill_;
      levels_below = std::max(levels_below, maximum_levels(next_univ));
    } else if (c.type_ == Fill::LATTICE) {
      Lattice& lat = *model::lattices[c.fill_];
      for (auto it = lat.begin(); it != lat.end(); ++it) {
        int32_t next_univ = *it;
        levels_below = std::max(levels_below, maximum_levels(next_univ));
      }
    }
  }

  ++levels_below;
  model::universe_level_counts[univ] = levels_below;
  return levels_below;
}

bool is_root_universe(int32_t univ_id)
{
  return model::universe_map[univ_id] == model::root_universe;
}

//==============================================================================

void free_memory_geometry()
{
  model::cells.clear();
  model::cell_map.clear();

  model::universes.clear();
  model::universe_map.clear();

  model::lattices.clear();
  model::lattice_map.clear();

  model::overlap_check_count.clear();
}

} // namespace openmc
