#ifndef OPENMC_RANDOM_RAY_TALLY_CONVERT_H
#define OPENMC_RANDOM_RAY_TALLY_CONVERT_H

namespace openmc {

// A mapping object that is used to map between a specific random ray 
// source region and an OpenMC native tally bin that it should score to
// every iteration.
struct TallyTask {
  int tally_idx;
  int filter_idx;
  int score_idx;
  int score_type;
  TallyTask(int tally_idx, int filter_idx, int score_idx, int score_type)
    : tally_idx(tally_idx), filter_idx(filter_idx), score_idx(score_idx), score_type(score_type)
  {}
};

namespace random_ray {

// 2D array representing values for all source regions x energy groups x tally tasks
extern std::vector<std::vector<TallyTask>> tally_task;

}

bool convert_source_regions_to_tallies();
void random_ray_tally();

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_TALLY_CONVERT_H
