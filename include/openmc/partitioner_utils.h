#ifndef OPENMC_PARTITIONER_UTILS_H
#define OPENMC_PARTITIONER_UTILS_H

#include "universe.h"
#include "position.h"
#include "timer.h"
#include "random_dist.h"

#include <omp.h>
#include <vector>
#include <set>

namespace openmc {

    using vec3 = Position;

    // This axis-aligned bounding box class (AABB) is designed to work easier with partitioners than openmc::BoundingBox
    struct AABB {
        AABB();
        AABB(const vec3& mi, const vec3& ma);

        void extend_max(const vec3& val);
        void extend_min(const vec3& val);

        void extend(const vec3& pos);
        void extend(const AABB& other_box);

        vec3 get_center() const;
        float surface_area() const;
        float volume() const;
        bool contains(const vec3& pos) const;

        void reset();

        bool operator==(const AABB& other) const;
        bool operator!=(const AABB& other) const;

        vec3 min;
        vec3 max;
    };

    template<class NodeT, size_t NUM_CHILDREN_PER_PARENT, size_t POOL_SIZE=16384>
    class NodeAllocator {
    public:
        NodeAllocator() : last_pool_next_index(0) {}

        ~NodeAllocator() {
            for(auto* ptr : pools) {
                delete[] ptr;
            }
        }

        NodeT* allocate() {
            if(last_pool_next_index == POOL_SIZE || pools.size() == 0) {
                pools.push_back(new NodeT[NUM_CHILDREN_PER_PARENT * POOL_SIZE]);
                last_pool_next_index = 0;
            }

            return &pools.back()[NUM_CHILDREN_PER_PARENT * last_pool_next_index++];
        }
    private:
        std::vector<NodeT*> pools;
        size_t last_pool_next_index;
    };

    struct CellPointUncompressed {
        vec3 pos;
        int cell;

        bool operator<(const CellPointUncompressed& other) const;
    };

    struct CellPoint {
        uint64_t data;

        int get_cell() const;
        int get_child_index(int depth) const;
        void compress_from(const CellPointUncompressed& uncomp, const AABB& bounds);

        bool operator<(const CellPoint& other) const;
    };

    class Bin {
    public:
        Bin();

        void insert(int cell);
        void insert_lockless(int cell);

        int size() const;
        int unique_size() const;

        void make_cells_unique();
        void sort_cells();

        float score(int iteration);

        const std::vector<int>& get_cells() const;

        void copy_untested_cells(const std::vector<int>& possible_cells, std::vector<int>& untested_cells);
    private:
        void lock_bin();
        void unlock_bin();

        omp_lock_t lock;
        std::vector<int> cells;
        int num_unique_cells;
        int prev_cell_count;
    };

    template<class NodeT>
    class ProbabilityBin {
    public:
        ProbabilityBin() : num_searched(0), num_found(0) {
            omp_init_lock(&found_lock);
        }

        ~ProbabilityBin() {
            omp_destroy_lock(&found_lock);

            for(auto& p : contained_nodes) {
                omp_destroy_lock(&p.second);
            }
        }

        float compute_score() {
            const float   REFINEMENT_BIN_SCORE_FOUND_MULT    = 1.0;
            const float   REFINEMENT_BIN_SCORE_SEARCHED_MULT = 1.0;

            float found_score = 
                  std::max(REFINEMENT_BIN_SCORE_FOUND_MULT * num_found, 1.0f) / 
                  std::max(REFINEMENT_BIN_SCORE_SEARCHED_MULT * num_searched, 1.0f);

            found_score *= found_score;

            float size_score = contained_nodes.size();
            return size_score;
        }

        size_t pick_node(uint64_t* seed) {
            return (size_t)uniform_distribution(0.0, (double)contained_nodes.size(), seed);
        }

        void update_common_cells() {
            found_cells.reserve(common_cells.size() + found_cells.size());
            
            for(int cell : common_cells) {
                found_cells.push_back(cell);
            }

            std::sort(found_cells.begin(), found_cells.end());

            common_cells.clear();

            int prev_cell = -1;
            for(int cell : found_cells) {
                if(cell != prev_cell) {
                    common_cells.push_back(cell);
                    prev_cell = cell;
                }
            }

            found_cells.clear();
        }

        void add_cell(int cell) {
            omp_set_lock(&found_lock);
            found_cells.push_back(cell);
            omp_unset_lock(&found_lock);
        }

        int num_searched, num_found;
        omp_lock_t found_lock;
        std::vector<int> common_cells;
        std::vector<int> found_cells;
        std::vector<std::pair<NodeT*, omp_lock_t>> contained_nodes;
    };

    template<typename StoreT, typename ToStoreT>
    void store_cell_point(StoreT& storage, const ToStoreT& to_store, const AABB& bounds);

    template<> inline
    void store_cell_point(CellPoint& storage, const CellPointUncompressed& to_store, const AABB& bounds) {
        storage.compress_from(to_store, bounds);
    }

    template<> inline
    void store_cell_point(CellPointUncompressed& storage, const CellPointUncompressed& to_store, const AABB& bounds) {
        storage = to_store;
    }

    template<typename T>
    std::vector<T> binned_point_search(const Universe& univ, const UniversePartitioner& fallback, const AABB& bounds)  {
        const int32_t BINNING_SEARCH_TOTAL_POINTS     = 8000000;
        const int32_t BINNING_SEARCH_TOTAL_ITERATIONS = 128;
        const int32_t BINNING_SEARCH_GRID_RES         = 32;

        float target_density = BINNING_SEARCH_TOTAL_POINTS / (BINNING_SEARCH_TOTAL_ITERATIONS * bounds.volume());

        std::vector<float> search_densities;
        for(int i = 0; i < BINNING_SEARCH_TOTAL_ITERATIONS; i++) { 
            search_densities.push_back(target_density);
        }

        // some console output vars
        int total_search_points = 0;
        for(float density : search_densities) {
            total_search_points += (int)(bounds.volume() * density);
        }

        std::vector<Bin> bin_grid(BINNING_SEARCH_GRID_RES * BINNING_SEARCH_GRID_RES * BINNING_SEARCH_GRID_RES); 
        std::vector<float> bin_cdf(BINNING_SEARCH_GRID_RES * BINNING_SEARCH_GRID_RES * BINNING_SEARCH_GRID_RES);

        vec3 bin_dim;
        for(int i = 0; i < 3; i++) {
            bin_dim[i] = (bounds.max[i] - bounds.min[i]) / BINNING_SEARCH_GRID_RES;
        }

        std::vector<T> points_in_bounds(total_search_points);
        int write_offset = 0;
        for(int iteration = 0; iteration < search_densities.size(); iteration++) {
            float density = search_densities[iteration];
            int num_search_points = (int)(bounds.volume() * density);

            float total_cdf = 0.0;
            for(int i = 0; i < bin_grid.size(); i++) {
                total_cdf += bin_grid[i].score(iteration);
                bin_cdf[i] = total_cdf;
            }

            for(float& cdf : bin_cdf) {
                cdf /= total_cdf;
            }

            #pragma omp parallel 
            {
                int tid = omp_get_thread_num();
                int tcount = omp_get_num_threads();

                uint64_t bin_seed = write_offset + tid;
                uint64_t pos_seed[3] = {
                    (uint64_t)(write_offset + tid + tcount), (uint64_t)(write_offset + tid + tcount * 2), (uint64_t)(write_offset + tid + tcount * 3)
                };

                std::vector<int> untested_cells;
                untested_cells.reserve(univ.cells_.size());

                #pragma omp for
                for(int i = 0; i < num_search_points; i++) {

                    int index = write_offset + i;
                    CellPointUncompressed point;

                    float bin_cdf_val = uniform_distribution(0.0, 1.0, &bin_seed);
                    auto cdf_iter = std::upper_bound(bin_cdf.begin(), bin_cdf.end(), bin_cdf_val);
                    int bin_idx = std::distance(bin_cdf.begin(), cdf_iter);
                    Bin* sample_bin = &bin_grid[bin_idx];

                    for(int j = 0; j < 3; j++) {
                        int idx = bin_idx % BINNING_SEARCH_GRID_RES;
                        bin_idx /= BINNING_SEARCH_GRID_RES;
                        float min = idx * bin_dim[j] + bounds.min[j];
                        float max = bin_dim[j] + min;

                        point.pos[j] = uniform_distribution(min, max, &pos_seed[j]);
                    }

                    Direction dummy_dir{1.0, 0.0, 0.0};

                    point.cell = univ.find_cell_for_point(sample_bin->get_cells(), point.pos);
                    if(point.cell == -1) {
                        const auto& possible_cells = fallback.get_cells(point.pos, dummy_dir);
                        sample_bin->copy_untested_cells(possible_cells, untested_cells);
                        point.cell = univ.find_cell_for_point(untested_cells, point.pos);

                        sample_bin->insert(point.cell);
                    } // else don't bother inserting

                    store_cell_point(points_in_bounds[index], point, bounds);
                }
            }
            
            #pragma omp parallel for
            for(auto& bin : bin_grid) {
                bin.sort_cells();
            }

            write_offset += num_search_points;
        }

        return points_in_bounds;  
    }

};

#endif // OPENMC_PARTITIONER_UTILS_H