#include "openmc/tallies/filter.h"

#include <fmt/core.h>

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/geometry_aux.h" // For distribcell_path
#include "openmc/lattice.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
Filter::DistribcellFilter_from_xml(pugi::xml_node node)
{
  auto cells = get_node_array<int32_t>(node, "bins");
  if (cells.size() != 1) {
    fatal_error("Only one cell can be specified per distribcell filter.");
  }

  // Find index in global cells vector corresponding to cell ID
  auto search = model::cell_map.find(cells[0]);
  if (search == model::cell_map.end()) {
    throw std::runtime_error{fmt::format(
      "Could not find cell {} specified on tally filter.", cell_)};
  }

  this->set_cell(search->second);
}

void
Filter::set_cell(int32_t cell)
{
  Expects(cell >= 0);
  Expects(cell < model::cells.size());
  cell_ = cell;
  n_bins_ = model::cells[cell].n_instances_;
}

void
Filter::DistribcellFilter_get_all_bins(const Particle& p, TallyEstimator estimator,
                                FilterMatch& match) const
{
  int offset = 0;
  auto distribcell_index = model::device_cells[cell_].distribcell_index_;
  for (int i = 0; i < p.n_coord_; i++) {
    auto& c {model::device_cells[p.coord_[i].cell]};
    if (c.type_ == Fill::UNIVERSE) {
      offset += c.offset_[distribcell_index];
    } else if (c.type_ == Fill::LATTICE) {
      auto& lat {model::device_lattices[p.coord_[i+1].lattice]};
      int i_xyz[3] {p.coord_[i+1].lattice_x,
                    p.coord_[i+1].lattice_y,
                    p.coord_[i+1].lattice_z};
      if (lat.are_valid_indices(i_xyz)) {
        offset += lat.offset(distribcell_index, i_xyz);
      }
    }
    if (cell_ == p.coord_[i].cell) {
      //match.bins_.push_back(offset);
      //match.weights_.push_back(1.0);
      match.bins_[match.bins_weights_length_] = offset;
      match.weights_[match.bins_weights_length_] = 1.0;
      match.bins_weights_length_++;
      return;
    }
  }
}

void
Filter::DistribcellFilter_to_statepoint(hid_t filter_group) const
{
  write_dataset(filter_group, "bins", model::cells[cell_].id_);
}

std::string
Filter::DistribcellFilter_text_label(int bin) const
{
  auto map = model::cells[cell_].distribcell_index_;
  auto path = distribcell_path(cell_, map, bin);
  return "Distributed Cell " + path;
}

} // namespace openmc
