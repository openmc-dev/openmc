#ifndef BINNING_H
#define BINNING_H

#include "universe.h"
#include "aabb.h"
#include <omp.h>
#include <vector>
#include <set>

namespace openmc {

    struct CellPointUncompressed {
        vec3 pos;
        int cell;
    };

    struct CellPoint {
        uint64_t data;

        int get_cell() const;
        int get_child_index(int depth) const;
        void compress_from(const CellPointUncompressed& uncomp, const AABB& bounds);
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

    std::vector<CellPoint> get_cell_points_binning(const Universe& univ, const UniversePartitioner& fallback, const AABB& bounds);
};

#endif