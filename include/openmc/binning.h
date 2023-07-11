#ifndef BINNING_H
#define BINNING_H

#include "universe.h"
#include "aabb.h"
#include "timer.h"
#include "random_dist.h"

#include <omp.h>
#include <vector>
#include <set>

namespace openmc {

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

        Timer t;
        t.start();
        // the purpose of CELL_CLEAR_FLAG is to mention when previous cells should be cleared
        // although we might want to keep all points, organizing the points into a tree baloons memory usage
        // hence we need to reduce points to prevent a crash
        // CELL_CLEAR_FLAG marks a location where we decide that the sampling method has learned enough and points before here can be discarded
        const float CELL_CLEAR_FLAG = -1.0;
        std::vector<float> K_SEARCH_DENSITY;// = {2.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0};

        float target_density = BINNING_SEARCH_TOTAL_POINTS / (BINNING_SEARCH_TOTAL_ITERATIONS * bounds.volume());

        K_SEARCH_DENSITY.clear();
        for(int i = 0; i < BINNING_SEARCH_TOTAL_ITERATIONS; i++) { 
            K_SEARCH_DENSITY.push_back(target_density);
        }


        // some console output vars
        int total_points_searched = 0;    
        omp_lock_t progess_display_lock;
        omp_init_lock(&progess_display_lock);

        int total_search_points = 0;
        for(float density : K_SEARCH_DENSITY) {
            if(density == CELL_CLEAR_FLAG)
                continue;

            total_search_points += (int)(bounds.volume() * density);
        }

        std::vector<Bin> bin_grid(BINNING_SEARCH_GRID_RES * BINNING_SEARCH_GRID_RES * BINNING_SEARCH_GRID_RES); 
        std::vector<float> bin_cdf(BINNING_SEARCH_GRID_RES * BINNING_SEARCH_GRID_RES * BINNING_SEARCH_GRID_RES);

        vec3 bin_dim;
        for(int i = 0; i < 3; i++) {
            bin_dim[i] = (bounds.max[i] - bounds.min[i]) / BINNING_SEARCH_GRID_RES;
        }

        std::cout << "Allocating " << total_search_points << " search points!" << std::endl;
        std::vector<T> points_in_bounds(total_search_points);
        int write_offset = 0;
        for(int iteration = 0; iteration < K_SEARCH_DENSITY.size(); iteration++) {
            float density = K_SEARCH_DENSITY[iteration];
            if(density == CELL_CLEAR_FLAG) {
                continue;
            }

            int num_search_points = (int)(bounds.volume() * density);

            float total_score = 0.0;
            for(int i = 0; i < bin_grid.size(); i++) {
                total_score += bin_grid[i].score(iteration);
                bin_cdf[i] = total_score;
            }

            for(int i = 0; i < bin_grid.size(); i++) {
                bin_cdf[i] /= total_score;
            }

            #pragma omp parallel 
            {
                int tid = omp_get_thread_num();
                int tcount = omp_get_num_threads();

                uint64_t bin_seed = write_offset + tid;
                uint64_t pos_seed[3] = {
                    (uint64_t)(tid + tcount), (uint64_t)(tid + tcount * 2), (uint64_t)(tid + tcount * 3)
                };

                std::vector<int> untested_cells;
                untested_cells.reserve(univ.cells_.size());

                #pragma omp for
                for(int i = 0; i < num_search_points; i++) {
                    #pragma omp atomic
                    total_points_searched++;

                    if(total_points_searched % 500000 == 0 && total_points_searched != 0) {
                        omp_set_lock(&progess_display_lock);
                        std::cout << "Currently searched " << total_points_searched << " out of " << total_search_points << " points (" << (100.0 * total_points_searched) / total_search_points << "%).\tETA is " << (float(total_search_points) / total_points_searched * t.elapsed()) - t.elapsed() << "\tmore seconds for this node.\n";
                        std::cout.flush();
                        omp_unset_lock(&progess_display_lock);
                    }

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

                    //points_in_bounds[index].compress_from(point, bounds);
                    store_cell_point(points_in_bounds[index], point, bounds);
                }
            }
            
            #pragma omp parallel for
            for(auto& bin : bin_grid) {
                bin.sort_cells();
            }

            write_offset += num_search_points;
        }

        int clear_pos = -1; 
        for(int i = K_SEARCH_DENSITY.size() - 1; i >= 0; i--) {
            if(K_SEARCH_DENSITY[i] == CELL_CLEAR_FLAG) {
                clear_pos = i;
                break;
            }
        }

        if(clear_pos != -1) {
            int total_discarded_search_points = 0;
            for(float density : K_SEARCH_DENSITY) {
                if(density == CELL_CLEAR_FLAG)
                    break;

                total_discarded_search_points += (int)(bounds.volume() * density);
            }

            points_in_bounds.erase(points_in_bounds.begin(), points_in_bounds.begin() + total_discarded_search_points);
            points_in_bounds.shrink_to_fit();
        }

        return points_in_bounds;  
    }
};

#endif