#include "openmc/binning.h"
#include "openmc/timer.h"
#include "openmc/random_dist.h"

namespace openmc {

    const int32_t CELL_POINT_COMPRESSION_MAX_DEPTH = 15;

    const int32_t BINNING_SEARCH_TOTAL_POINTS     = 8000000;
    const int32_t BINNING_SEARCH_TOTAL_ITERATIONS = 128;
    const int32_t BINNING_SEARCH_GRID_RES         = 32;

    int CellPoint::get_cell() const {
        uint64_t index = (data >> (3 * CELL_POINT_COMPRESSION_MAX_DEPTH));
        return (int)index;
    }

    int CellPoint::get_child_index(int depth) const {
        uint64_t index = ((data >> (3 * depth)) & 0b111);
        return (int)index;
    }

    void CellPoint::compress_from(const CellPointUncompressed& uncomp, const AABB& bounds) {
        data = 0;

        vec3 node_dim;
        for(int i = 0; i < 3; i++) {
            node_dim[i] = (bounds.max[i] - bounds.min[i]) * 0.5;
        }

        vec3 center = bounds.get_center();

        for(int i = 0; i < CELL_POINT_COMPRESSION_MAX_DEPTH; i++) {
            // halve the node dim
            uint64_t idx = 0;
            for(int i = 0; i < 3; i++) {
                node_dim[i] *= 0.5;
                bool less = (uncomp.pos[i] < center[i]);
                center[i] += node_dim[i] * (less ? -1 : 1);
                idx = 2 * idx + int(less);
            }

            idx = (idx << (3 * i));
            data |= idx;
        }

        uint64_t cell_comp = uncomp.cell;
        cell_comp = (cell_comp << (3 * CELL_POINT_COMPRESSION_MAX_DEPTH));
        data |= cell_comp;
    }

    Bin::Bin() : num_unique_cells(0), prev_cell_count(0) {
        omp_init_lock(&lock);
        cells.reserve(512);
    }

    void Bin::insert(int cell) {
        lock_bin();
        insert_lockless(cell);
        unlock_bin();
    }

    void Bin::insert_lockless(int cell) {
        cells.push_back(cell);
        if(cells.size() == cells.capacity()) {
            int cur_size = cells.size();

            make_cells_unique();

            int size_saving = cur_size - cells.size();
            if(size_saving < (int)(0.25 * cur_size)) {
                cells.reserve(2 * cells.capacity());
            }
        }            
    }

    void Bin::sort_cells() {
        std::sort(cells.begin(), cells.end());
    }

    void Bin::make_cells_unique() {
        sort_cells();
        int i = 0, j = 1;
        while(j < cells.size()) {
            if(cells[j] != cells[i]) {
                i++;
                cells[i] = cells[j];
            }
            j++;
        }
        cells.resize(i + 1);
        num_unique_cells = cells.size();
    }

    int Bin::size() const {
        return std::max(cells.size(), (size_t)1);
    }

    int Bin::unique_size() const {
        return num_unique_cells;
    }

    float Bin::score(int iteration) {
        float val;

        if(iteration == 0) {
            val = 1.0;
        } else if(iteration == 1) {
            val = (float)size();
        } else {
            const float alpha = 0.1;

            float size_score = size();
            float cell_increase_score = float(cells.size()) / prev_cell_count;

            val = alpha * cell_increase_score + (1 - alpha) * size_score;
        }

        // do stuff here to update for the next scoring cycle
        prev_cell_count = size();

        return val;
    }

    const std::vector<int>& Bin::get_cells() const {
        return cells;
    }

    void Bin::copy_untested_cells(const std::vector<int>& possible_cells, std::vector<int>& untested_cells) {
        untested_cells.clear();

        lock_bin();
        for(int cell : possible_cells) {
            if(!std::binary_search(cells.begin(), cells.begin() + num_unique_cells, cell)) {
                untested_cells.push_back(cell);
            }
        }
        unlock_bin();

    }
    void Bin::lock_bin() {
        omp_set_lock(&lock);
    }

    void Bin::unlock_bin() {
        omp_unset_lock(&lock);
    }
    
    std::vector<CellPoint> get_cell_points_binning(const Universe& univ, const UniversePartitioner& fallback, const AABB& bounds) {
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
        std::vector<CellPoint> points_in_bounds(total_search_points);
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

                    points_in_bounds[index].compress_from(point, bounds);
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