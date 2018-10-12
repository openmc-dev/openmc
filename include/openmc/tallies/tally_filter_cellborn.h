#ifndef OPENMC_TALLY_FILTER_CELLBORN_H
#define OPENMC_TALLY_FILTER_CELLBORN_H

#include "openmc/tallies/tally_filter_cell.h"

namespace openmc {

class CellbornFilter : public CellFilter
{
public:
  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
  {
    auto search = map_.find(p->cell_born);
    if (search != map_.end()) {
      // TODO: off-by-one
      match.bins.push_back(search->second + 1);
      match.weights.push_back(1);
    }
  }

  virtual void
  to_statepoint(hid_t filter_group) const override
  {
    write_dataset(filter_group, "type", "cellborn");
    write_dataset(filter_group, "n_bins", n_bins_);
    std::vector<int32_t> cell_ids;
    for (auto c : cells_) cell_ids.push_back(cells[c]->id_);
    write_dataset(filter_group, "bins", cell_ids);
  }

  virtual std::string
  text_label(int bin) const override
  {
    return "Birth Cell " + std::to_string(cells[cells_[bin-1]]->id_);
  }
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_CELLBORN_H
