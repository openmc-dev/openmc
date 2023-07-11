#include "openmc/binning.h"


namespace openmc {

    const int32_t CELL_POINT_COMPRESSION_MAX_DEPTH = 15;

    bool CellPointUncompressed::operator<(const CellPointUncompressed& other) const {
        return (cell < other.cell);
    }

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

    bool CellPoint::operator<(const CellPoint& other) const {
        return (get_cell() < other.get_cell());
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

};