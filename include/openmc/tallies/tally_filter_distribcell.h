#ifndef OPENMC_TALLY_FILTER_DISTRIBCELL_H
#define OPENMC_TALLY_FILTER_DISTRIBCELL_H

#include <string>

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/geometry_aux.h" // For distribcell_path
#include "openmc/lattice.h"
#include "openmc/tallies/tally_filter.h"


namespace openmc {

//==============================================================================
//! Specifies which distributed geometric cells tally events reside in.
//==============================================================================

class DistribcellFilter : public TallyFilter
{
public:
  std::string type() const override {return "distribcell";}

  ~DistribcellFilter() = default;

  void
  from_xml(pugi::xml_node node) override
  {
    auto cells = get_node_array<int32_t>(node, "bins");
    if (cells.size() != 1) {
      fatal_error("Only one cell can be specified per distribcell filter.");
    }
    cell_ = cells[0];
  }

  void
  initialize() override
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
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
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
  to_statepoint(hid_t filter_group) const override
  {
    TallyFilter::to_statepoint(filter_group);
    write_dataset(filter_group, "bins", cells[cell_]->id_);
  }

  std::string
  text_label(int bin) const override
  {
    auto map = cells[cell_]->distribcell_index_;
    auto path = distribcell_path(cell_, map, bin-1);

    return "Distributed Cell " + path;
  }

  int32_t cell_;
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_DISTRIBCELL_H
