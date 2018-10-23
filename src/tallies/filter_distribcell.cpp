#include "openmc/tallies/filter_distribcell.h"

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/geometry_aux.h" // For distribcell_path
#include "openmc/lattice.h"

namespace openmc {

void
DistribcellFilter::from_xml(pugi::xml_node node)
{
  auto cells = get_node_array<int32_t>(node, "bins");
  if (cells.size() != 1) {
    fatal_error("Only one cell can be specified per distribcell filter.");
  }
  cell_ = cells[0];
}

void
DistribcellFilter::initialize()
{
  auto search = cell_map.find(cell_);
  if (search != cell_map.end()) {
    cell_ = search->second;
    n_bins_ = cells[cell_]->n_instances_;
  } else {
    std::stringstream err_msg;
    err_msg << "Could not find cell " << cell_
            << " specified on tally filter.";
    fatal_error(err_msg);
  }
}

void
DistribcellFilter::get_all_bins(Particle* p, int estimator, FilterMatch& match)
const
{
  int offset = 0;
  auto distribcell_index = cells[cell_]->distribcell_index_;
  for (int i = 0; i < p->n_coord; i++) {
    auto& c {*cells[p->coord[i].cell]};
    if (c.type_ == FILL_UNIVERSE) {
      offset += c.offset_[distribcell_index];
    } else if (c.type_ == FILL_LATTICE) {
      auto& lat {*lattices[p->coord[i+1].lattice-1]};
      int i_xyz[3] {p->coord[i+1].lattice_x,
                    p->coord[i+1].lattice_y,
                    p->coord[i+1].lattice_z};
      if (lat.are_valid_indices(i_xyz)) {
        offset += lat.offset(distribcell_index, i_xyz);
      }
    }
    if (cell_ == p->coord[i].cell) {
      match.bins_.push_back(offset + 1);
      match.weights_.push_back(1);
      return;
    }
  }
}

void
DistribcellFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", cells[cell_]->id_);
}

std::string
DistribcellFilter::text_label(int bin) const
{
  auto map = cells[cell_]->distribcell_index_;
  auto path = distribcell_path(cell_, map, bin-1);
  return "Distributed Cell " + path;
}

} // namespace openmc
