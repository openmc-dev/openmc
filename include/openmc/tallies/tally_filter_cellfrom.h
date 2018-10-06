#ifndef OPENMC_TALLY_FILTER_CELLFROM_H
#define OPENMC_TALLY_FILTER_CELLFROM_H

#include "openmc/tallies/tally_filter_cell.h"


namespace openmc {

class CellFromFilter : public CellFilter
{
public:
  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
  {
    for (int i = 0; i < p->last_n_coord; i++) {
      auto search = map_.find(p->last_cell[i]);
      if (search != map_.end()) {
        // TODO: off-by-one
        match.bins.push_back(search->second + 1);
        match.weights.push_back(1);
      }
    }
  }

  virtual void
  to_statepoint(hid_t filter_group) const override
  {
    write_dataset(filter_group, "type", "cellfrom");
    write_dataset(filter_group, "n_bins", cells_.size());
    std::vector<int32_t> cell_ids;
    for (auto c : cells_) cell_ids.push_back(cells[c]->id_);
    write_dataset(filter_group, "bins", cell_ids);
  }

  virtual std::string
  text_label(int bin) const override
  {
    return "Cell from " + std::to_string(cells[cells_[bin-1]]->id_);
  }
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_CELLFROM_H
