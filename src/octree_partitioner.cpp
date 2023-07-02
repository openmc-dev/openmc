#include "openmc/octree_partitioner.h"
#include "openmc/aabb.h"
#include "openmc/random_dist.h"
#include "openmc/timer.h"

#include <stack>
#include <queue>
#include <thread>
#include <algorithm>
#include <assert.h>
#include <omp.h>

#include <fstream>

#define OCTREE_DEBUG_INFO

namespace openmc {

struct CellPoint {
    vec3 pos;
    int cell;
};

int octree_use_count = 0;
int fallback_use_count = 0;

#ifdef OCTREE_DEBUG_INFO
omp_lock_t ocm; // short for octree console mutex
bool ocm_init = false; // unlikely that it will be double init

// shorten some stuff
#define LOCK_OCM() omp_set_lock(&ocm)
#define UNLOCK_OCM() omp_unset_lock(&ocm)
#endif

OctreeNode::OctreeNode() : children(nullptr), depth(0) {}

bool OctreeNode::is_leaf() const {
    return (children == nullptr);
}

int OctreeNode::get_containing_child_index(const vec3& r) const {
    int idx = 4 * int(r.x < center.x) + 2 * int(r.y < center.y) + int(r.z < center.z);
    return idx;
}

OctreeNode& OctreeNode::get_containing_child(const vec3& r) const {
    return children[get_containing_child_index(r)];
}

std::vector<AABB> OctreeNode::subdivide(const AABB& parent) {
    std::vector<AABB> resultant_boxes;
    resultant_boxes.push_back(parent);
    for(int i = 0; i < 3; i++) {
        std::vector<AABB> temp_box_buffer;

        for(const auto& box : resultant_boxes) {
            // split on i-th axis
            float midpoint = box.get_center()[i];

            int j = ((i + 1) % 3);
            int k = ((i + 2) % 3);

            AABB splitted_boxes[2];
            vec3 extension_point[2];

            extension_point[0][i] = midpoint;
            extension_point[0][j] = box.min[j];
            extension_point[0][k] = box.min[k];

            extension_point[1][i] = midpoint;
            extension_point[1][j] = box.max[j];
            extension_point[1][k] = box.max[k];

            splitted_boxes[0].extend(box.max);
            splitted_boxes[0].extend(extension_point[0]);

            splitted_boxes[1].extend(box.min);
            splitted_boxes[1].extend(extension_point[1]);

            temp_box_buffer.push_back(splitted_boxes[0]);
            temp_box_buffer.push_back(splitted_boxes[1]);
        }

        // move the results to the next splitting stage
        resultant_boxes = temp_box_buffer;
    }

    for(int i = 0; i < 8; i++) {
        children[i].center = resultant_boxes[i].get_center();
    }

    return resultant_boxes;
}
 
void find_cells_in_box(const Universe& univ, const OctreeNode& parent, OctreeNode& child, const AABB& box) {
    const double NUM_CELL_SEARCH_POINT_DENSITY = 4.0;
    const int num_search_points = (int)(box.volume() * NUM_CELL_SEARCH_POINT_DENSITY);

    // get some numbers that are unique across all invocations
    uint64_t octree_prng_seed[3] = {(uint64_t)&child, (uint64_t)child.id, (uint64_t)((child.id * 3244230925 + 432534) % 436786)};

    #ifdef OCTREE_DEBUG_INFO
    const int POINTS_LOGGING_THRESHOLD = 16384; 
    const int POINTS_LOGGING_INCREMENT = 5000;
    if(num_search_points > POINTS_LOGGING_THRESHOLD) {
        LOCK_OCM();
        std::cout << "Note: searching node " << child.id << " for " << num_search_points << " points.\n";
        UNLOCK_OCM();
    }

    Timer t;
    t.start();
    #endif

    std::vector<bool> skip_cell(parent.cells.size());
    child.cells.reserve(parent.cells.size());
    for(int i = 0; i < num_search_points; i++) {
        #ifdef OCTREE_DEBUG_INFO
        if(num_search_points > POINTS_LOGGING_THRESHOLD && i % POINTS_LOGGING_INCREMENT == 0 && i != 0) {
            LOCK_OCM();
            std::cout << "Node " << child.id << ":\tcurrently searched " << i << " out of " << num_search_points << " points (" << (100.0 * i) / num_search_points 
                      << "%) and found " << child.cells.size() << " unique cells. ";
            std::cout << "ETA is " << ((float)num_search_points / i * t.elapsed()) - t.elapsed() << " more seconds for this node.\n";
            std::cout.flush(); // you may want to comment this out
            UNLOCK_OCM();
        }
        #endif

        vec3 rand_pos;

        // gen random point
        rand_pos.x = uniform_distribution(box.min.x, box.max.x, &octree_prng_seed[0]);
        rand_pos.y = uniform_distribution(box.min.y, box.max.y, &octree_prng_seed[1]);
        rand_pos.z = uniform_distribution(box.min.z, box.max.z, &octree_prng_seed[2]);

        univ.find_cell_in_list(parent.cells, child.cells, skip_cell, rand_pos);
    }
} 

OctreeNodeAllocator::OctreeNodeAllocator() : last_pool_next_index(0) {}

OctreeNodeAllocator::~OctreeNodeAllocator() {
    for(auto* ptr : pools) {
        delete[] ptr;
    }
}

const int OCTREE_NODE_ALLOCATOR_POOL_SIZE = 1024; // number of 8 children packs
OctreeNode* OctreeNodeAllocator::allocate() {
    if(last_pool_next_index == OCTREE_NODE_ALLOCATOR_POOL_SIZE || pools.size() == 0) {
        pools.push_back(new OctreeNode[8 * OCTREE_NODE_ALLOCATOR_POOL_SIZE]);
        last_pool_next_index = 0;
    }

    auto ptr = &pools.back()[8 * last_pool_next_index];
    last_pool_next_index++;
    return ptr;
}

Timer t;
std::vector<CellPoint> get_cell_points_uniform_dist(const Universe& univ, const UniversePartitioner& fallback, const AABB& bounds) {
    const float NUM_CELL_SEARCH_POINT_DENSITY = 32.0f;
    const int num_search_points = (int)(bounds.volume() * NUM_CELL_SEARCH_POINT_DENSITY);
    std::vector<CellPoint> points_in_bounds(num_search_points);

    int total_points_searched = 0;
        
    omp_lock_t progess_display_lock;
    omp_init_lock(&progess_display_lock);

    #pragma omp parallel for
    for(int i = 0; i < num_search_points; i++) {
        #pragma omp atomic
        total_points_searched++;

        if(total_points_searched % 500000 == 0 && total_points_searched != 0) {
            omp_set_lock(&progess_display_lock);
            std::cout << "Currently searched " << total_points_searched << " out of " << num_search_points << " points (" << (100.0 * total_points_searched) / num_search_points << "%).\tETA is " << (float(num_search_points) / total_points_searched * t.elapsed()) - t.elapsed() << "\tmore seconds for this node.\n";
            std::cout.flush();
            omp_unset_lock(&progess_display_lock);
        }

        auto& point = points_in_bounds[i];

        uint64_t seed[3] = {
            (uint64_t)i, (uint64_t)(num_search_points + i), (uint64_t)(2 * num_search_points + i)
        };

        for(int j = 0; j < 3; j++) {
            point.pos[j] = uniform_distribution(bounds.min[j], bounds.max[j], &seed[j]);
        }

        Direction dummy_dir{1.0, 0.0, 0.0};

        const auto& possible_cells = fallback.get_cells(point.pos, dummy_dir);
        point.cell = univ.find_cell_for_point(possible_cells, point.pos);
    }

    return points_in_bounds;
}



std::vector<CellPoint> get_cell_points_binning(const Universe& univ, const UniversePartitioner& fallback, const AABB& bounds) {
    // the purpose of CELL_CLEAR_FLAG is to mention when previous cells should be cleared
    // although we might want to keep all points, organizing the points into a tree baloons memory usage
    // hence we need to reduce points to prevent a crash
    // CELL_CLEAR_FLAG marks a location where we decide that the sampling method has learned enough and points before here can be discarded
    const float CELL_CLEAR_FLAG = -1.0;
    const std::vector<float> NUM_CELL_SEARCH_POINT_DENSITY = {
        2.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0
    };

    /*
    NUM_CELL_SEARCH_POINT_DENSITY.clear();
    for(int i = 0; i < 4; i++) {
        NUM_CELL_SEARCH_POINT_DENSITY.push_back(8.0);
    }
    */

    const int BIN_GRID_RES = 32;

    // some console output vars
    int total_points_searched = 0;    
    omp_lock_t progess_display_lock;
    omp_init_lock(&progess_display_lock);

    int total_search_points = 0;
    for(float density : NUM_CELL_SEARCH_POINT_DENSITY) {
        if(density == CELL_CLEAR_FLAG)
            continue;

        total_search_points += (int)(bounds.volume() * density);
    }

    #define IN_PLACE_UNIQUIFIYING // in place unique-ifiying 
    const int RESERVATION_SIZE = 512;
    class Bin {
    public:
        Bin() : num_unique_cells(0), prev_cell_count(0) {
            omp_init_lock(&lock);
            cells.reserve(RESERVATION_SIZE);
        }

        void insert(const CellPoint& p) {
            lock_bin();
            insert_lockless(p);
            unlock_bin();
        }

        void insert_lockless(const CellPoint& p) {
            cells.push_back(p.cell);
            if(cells.size() == cells.capacity()) {
                int cur_size = cells.size();

                make_cells_unique();

                int size_saving = cur_size - cells.size();
                if(size_saving < (int)(0.25 * cur_size)) {
                    cells.reserve(2 * cells.capacity());
                }
            }            
        }

        void sort_cells() {
            std::sort(cells.begin(), cells.end());
        }

        void make_cells_unique() {
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

        int size() const {
            return std::max(cells.size(), (size_t)1);
        }

        int unique_size() const {
            return num_unique_cells;
        }

        float score(int iteration) {
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

        const std::vector<int>& get_cells() const {
            return cells;
        }

        void copy_untested_cells(const std::vector<int>& possible_cells, std::vector<int>& untested_cells) {
            untested_cells.clear();

            lock_bin();
            for(int cell : possible_cells) {
                if(!std::binary_search(cells.begin(), cells.begin() + num_unique_cells, cell)) {
                    untested_cells.push_back(cell);
                }
            }
            unlock_bin();

        }
    private:
        void lock_bin() {
            omp_set_lock(&lock);
        }

        void unlock_bin() {
            omp_unset_lock(&lock);
        }

        omp_lock_t lock;
        std::vector<int> cells;
        int num_unique_cells;
        int prev_cell_count;
 
    };

    std::vector<Bin> bin_grid(BIN_GRID_RES * BIN_GRID_RES * BIN_GRID_RES); 
    std::vector<float> bin_cdf(BIN_GRID_RES * BIN_GRID_RES * BIN_GRID_RES);

    vec3 bin_dim;
    for(int i = 0; i < 3; i++) {
        bin_dim[i] = (bounds.max[i] - bounds.min[i]) / BIN_GRID_RES;
    }

    std::vector<CellPoint> points_in_bounds(total_search_points);
    int write_offset = 0;
    for(int iteration = 0; iteration < NUM_CELL_SEARCH_POINT_DENSITY.size(); iteration++) {
        float density = NUM_CELL_SEARCH_POINT_DENSITY[iteration];
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
                auto& point = points_in_bounds[index];

                float bin_cdf_val = uniform_distribution(0.0, 1.0, &bin_seed);
                auto cdf_iter = std::upper_bound(bin_cdf.begin(), bin_cdf.end(), bin_cdf_val);
                int bin_idx = std::distance(bin_cdf.begin(), cdf_iter);
                Bin* sample_bin = &bin_grid[bin_idx];

                for(int j = 0; j < 3; j++) {
                    int idx = bin_idx % BIN_GRID_RES;
                    bin_idx /= BIN_GRID_RES;
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

                    sample_bin->insert(point);
                } // else don't bother inserting
            }
        }
        
        #pragma omp parallel for
        for(auto& bin : bin_grid) {
            bin.sort_cells();
        }

        write_offset += num_search_points;
    }

    int clear_pos = -1; 
    for(int i = NUM_CELL_SEARCH_POINT_DENSITY.size() - 1; i >= 0; i--) {
        if(NUM_CELL_SEARCH_POINT_DENSITY[i] == CELL_CLEAR_FLAG) {
            clear_pos = i;
            break;
        }
    }

    if(clear_pos != -1) {
        int total_discarded_search_points = 0;
        for(float density : NUM_CELL_SEARCH_POINT_DENSITY) {
            if(density == CELL_CLEAR_FLAG)
                break;

            total_discarded_search_points += (int)(bounds.volume() * density);
        }

        points_in_bounds.erase(points_in_bounds.begin(), points_in_bounds.begin() + total_discarded_search_points);
        points_in_bounds.shrink_to_fit();
    }

    return points_in_bounds;
}

OctreePartitioner::OctreePartitioner(const Universe& univ, int target_cells_per_node, int max_depth, const std::string& file_path) : fallback(univ) {
    if(!ocm_init) {
        ocm_init = true;
        omp_init_lock(&ocm);
    }

    t.start();
    std::cout << "====================OCTREE CONSTRUCTION===================\n";
    int max_tree_size = (2 << (3 * max_depth));
    std::cout << "Maximum tree size is " << max_tree_size << '\n';
    // octrees don't have to have cube nodes/leaves, so we don't have to convert the source AABB into a square
    auto urootbox = univ.bounding_box();
    AABB rootbox(vec3(urootbox.xmin, urootbox.ymin, urootbox.zmin), vec3(urootbox.xmax, urootbox.ymax, urootbox.zmax));

    const float half_side_length = 50.0;
    rootbox.min = vec3(-half_side_length, -half_side_length, -half_side_length);
    rootbox.max = vec3( half_side_length,  half_side_length,  half_side_length);

    bounds = rootbox;
    root.center = rootbox.get_center();
    //root.cells = univ.cells_;
    root.depth = 0;
    root.id = 0;

    if(univ.cells_.size() <= target_cells_per_node) {
        std::cout << "Universe has only " << root.cells.size() << " cells, which is below the target cells per node, " 
                  << target_cells_per_node << ". Octree will have depth of only 1.\n";
        root.cells = univ.cells_;
        return;
    }



    std::vector<CellPoint> points_in_bounds = get_cell_points_binning(univ, fallback, bounds);

    std::cout << "Done searching for points! Beginning organization of points in octree...\n";
    std::cout.flush();

    struct OctreeConstructionTask {
        OctreeNode* node;
        AABB box;
        std::vector<CellPoint> points;

        OctreeConstructionTask() = default;
        OctreeConstructionTask(OctreeNode* n, const AABB& b, const std::vector<CellPoint>& p) : node(n), box(b), points(p) {}
    };

    // kind of a crude way to sort but whatever
    Timer sorting_timer;
    sorting_timer.start();
   
    struct PointComp {
        inline bool operator() (const CellPoint& lhs, const CellPoint& rhs) {
            return (lhs.cell < rhs.cell);
        }
    }; 
    std::sort(points_in_bounds.begin(), points_in_bounds.end(), PointComp());

    std::cout << "Total sorting time: " << sorting_timer.elapsed() << '\n';

    double popping_time = 0.0, init_time = 0.0, alloc_time = 0.0, node_selec_time = 0.0, post_process_time = 0.0;

    std::stack<OctreeConstructionTask> oct_unproc;
    oct_unproc.emplace(&root, bounds, points_in_bounds);
    Timer total_building_timer;
    total_building_timer.start();
    while(!oct_unproc.empty()) {
        Timer t;

        t.start();
        auto cur_task = std::move(oct_unproc.top());
        oct_unproc.pop();
        popping_time += t.elapsed();
        t.reset();

        // subdivide
        t.start();

        cur_task.node->children = node_alloc.allocate();

        auto boxes = cur_task.node->subdivide(cur_task.box);
        init_time += t.elapsed();
        t.reset();

        // allocation
        t.start();
        OctreeConstructionTask child_tasks[8];
        for(int i = 0; i < 8; i++) {
            child_tasks[i].box = boxes[i];
            child_tasks[i].node = &cur_task.node->children[i];
            child_tasks[i].points.reserve(cur_task.points.size());
        }
        alloc_time += t.elapsed();
        t.reset();

        // sort points
        t.start();
        for(const auto& point : cur_task.points) {
            child_tasks[cur_task.node->get_containing_child_index(point.pos)].points.push_back(point);
        }
        node_selec_time += t.elapsed();
        t.reset();

        // post processing (make nodes leaves or push on construction stack)
        t.start();
        for(int i = 0; i < 8; i++) {
            // we don't have to sort since all points were already sorted before processing loop
            int prev_cell = -1;
            int num_unique_cells = 0;
            for(const auto& p : child_tasks[i].points) {
                if(p.cell != prev_cell) {
                    num_unique_cells++;
                    prev_cell = p.cell;
                }
            }


            if(num_unique_cells > target_cells_per_node) {
                oct_unproc.push(std::move(child_tasks[i]));
            } else {
                // now, make the points unique
                
                child_tasks[i].node->cells.resize(num_unique_cells);

                prev_cell = -1;
                for(const auto& p : child_tasks[i].points) {
                    if(p.cell != prev_cell) {
                        child_tasks[i].node->cells.push_back(p.cell);
                        prev_cell = p.cell;
                    }
                }               
            }
        }
        post_process_time += t.elapsed();
        t.reset();
    }
    total_building_timer.stop();

    std::cout << "Done placing points into octree! Performance statistics:\n";
    std::cout << "\tTask popping:\t" << popping_time << '\n';
    std::cout << "\tChild init:\t" << init_time << '\n';
    std::cout << "\tAllocation:\t" << alloc_time << '\n';
    std::cout << "\tNode selec:\t" << node_selec_time << '\n';
    std::cout << "\tPost proc:\t" << post_process_time << '\n';
    std::cout << "\tExpec tot:\t" << popping_time + init_time + alloc_time + node_selec_time + post_process_time << '\n';
    std::cout << "\tTotal time:\t" << total_building_timer.elapsed() << '\n';


    t.stop();
    std::cout << "Construction took " << t.elapsed() << " seconds.\n";

    // now, we write the octree to disk if a path is specified
    if(file_path.size() != 0) {
        write_to_file(file_path);
    }

}  

OctreePartitioner::OctreePartitioner(const Universe& univ, const std::string& file_path) : fallback(univ) {
    read_from_file(file_path);
}

OctreePartitioner::~OctreePartitioner() {
    std::cout << "Successful octree uses:\t" << octree_use_count << '\n';
    std::cout << "Octree fallback uses:\t" << fallback_use_count << '\n';
    std::cout << "Again, construction time:\t" << t.elapsed() << '\n';
}

//#define FALLBACK_OPTIMIZATIONS

const vector<int32_t>& OctreePartitioner::get_cells(Position r, Direction u) const {
    octree_use_count++;

    if(!bounds.contains(r)) {
        // discount this get_cells call from the stats to only focus on points within the octree
        fallback_use_count--; 
        octree_use_count--;
        return get_cells_fallback(r, u);
    }


    const OctreeNode* current = &root;

    while(!current->is_leaf()) {
        current = &current->get_containing_child(r);
    }

    if(current->cells.empty()) {
        return get_cells_fallback(r, u);
    }

    return current->cells;
}

const vector<int32_t>& OctreePartitioner::get_cells_fallback(Position r, Direction u) const {
    fallback_use_count++;

    return fallback.get_cells(r, u);
}

struct OctreeNodeSerialized {
    int id;

    vec3 center;
    bool is_leaf;

    // parent data
    int first_child_index;

    // leaf data
    int contained_cells_index;
    int num_contained_cells;
};

void OctreePartitioner::write_to_file(const std::string& file_path){
    std::fstream octree_file(file_path, std::ios::out | std::ios::binary);

    std::vector<OctreeNodeSerialized> node_data;
    std::vector<int> cell_data;

    // serialize
    std::queue<const OctreeNode*> unproc_nodes;
    unproc_nodes.push(&root);

    while(!unproc_nodes.empty()) {
        auto& cur = *unproc_nodes.front();
        unproc_nodes.pop();

        OctreeNodeSerialized ser;

        ser.id = cur.id;
        ser.center = cur.center;

        if(ser.is_leaf = cur.is_leaf()) {
            ser.num_contained_cells = cur.cells.size();
            ser.contained_cells_index = cell_data.size();

            for(int cell : cur.cells) {
                cell_data.push_back(cell);
            }
        } else {
            ser.first_child_index = 1 + node_data.size() + unproc_nodes.size();
            for(int i = 0; i < 8; i++) {
                unproc_nodes.push(&cur.children[i]);
            }
        }
        
        node_data.push_back(ser);
    }

    // write to disk
    #define WR_BIN(x) octree_file.write((char*)&x, sizeof(x))

    WR_BIN(num_nodes);
    for(const auto& node : node_data) {
        WR_BIN(node);
    }    

    int num_cells = cell_data.size();
    WR_BIN(num_cells);
    for(auto cell : cell_data) {
        WR_BIN(cell);
    }
}

void OctreePartitioner::read_from_file(const std::string& file_path) {
    std::fstream octree_file(file_path, std::ios::in | std::ios::binary);

    std::vector<OctreeNodeSerialized> node_data;
    std::vector<int> cell_data;

    #define RD_BIN(x) octree_file.read((char*)&x, sizeof(x))

    RD_BIN(num_nodes);

    node_data.resize(num_nodes);
    for(auto& node : node_data) {
        RD_BIN(node);
    }

    int num_cells;
    RD_BIN(num_cells);
    cell_data.resize(num_cells);
    for(auto& cell : cell_data) {
        RD_BIN(cell);
    }

    std::stack<std::pair<OctreeNode*, int>> unproc_nodes;
    unproc_nodes.emplace(&root, 0);

    while(!unproc_nodes.empty()) {
        auto& p = unproc_nodes.top();
        unproc_nodes.pop();

        OctreeNode& cur = *p.first;
        int index = p.second;

        OctreeNodeSerialized& ser = node_data[index];

        cur.id = ser.id;
        cur.center = ser.center;

        if(ser.is_leaf) {
            cur.cells.reserve(ser.num_contained_cells);
            for(int i = ser.contained_cells_index; i < ser.contained_cells_index + ser.num_contained_cells; i++) {
                cur.cells.push_back(cell_data[i]);
            }
        } else {
            cur.children = new OctreeNode[8];

            for(int i = 0; i < 8; i++) {
                unproc_nodes.emplace(&cur.children[i], ser.first_child_index + i);
            }
        }
    }
}

}

/*
parallel splitting algorithm:

Successful octree uses:	120035
Octree fallback uses:	42311
Again, construction time:	126.471
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 23.7694 seconds.
 Universe::find_cell took 5.79047 seconds in total. (Approximately 24.361% of total execution time)


Successful octree uses:	120181
Octree fallback uses:	72779
Again, construction time:	59.4075
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 26.0166 seconds.
 Universe::find_cell took 7.97749 seconds in total. (Approximately 30.663% of total execution time)


Successful octree uses:	120036
Octree fallback uses:	55888
Again, construction time:	95.169
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 24.5633 seconds.
 Universe::find_cell took 6.71337 seconds in total. (Approximately 27.3309% of total execution time)


Successful octree uses:	120072
Octree fallback uses:	42043
Again, construction time:	126.448
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 24.1136 seconds.
 Universe::find_cell took 5.78582 seconds in total. (Approximately 23.994% of total execution time)


alpha 0.1, 2-2-4-8-16
Successful octree uses:	120078
Octree fallback uses:	42061
Again, construction time:	123.161
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 23.8965 seconds.
 Universe::find_cell took 5.80047 seconds in total. (Approximately 24.2733% of total execution time)


alpha 0.4
Successful octree uses:	119959
Octree fallback uses:	45403
Again, construction time:	116.55
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 24.2023 seconds.
 Universe::find_cell took 6.06049 seconds in total. (Approximately 25.041% of total execution time)


without in-place op
Successful octree uses:	120009
Octree fallback uses:	44386
Again, construction time:	120.182
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 24.3735 seconds.
 Universe::find_cell took 6.08325 seconds in total. (Approximately 24.9585% of total execution time)


after the fix (4-12-32)
Successful octree uses:	120090
Octree fallback uses:	44398
Again, construction time:	119.082
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 24.2427 seconds.
 Universe::find_cell took 5.94023 seconds in total. (Approximately 24.5032% of total execution time)


with binning
Successful octree uses:	120001
Octree fallback uses:	72621
Again, construction time:	77.2118
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 25.9852 seconds.
 Universe::find_cell took 7.89101 seconds in total. (Approximately 30.3673% of total execution time)


with no bining, constant dens 128
Successful octree uses:	120034
Octree fallback uses:	46379
Again, construction time:	233.664
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 24.2725 seconds.
 Universe::find_cell took 6.091 seconds in total. (Approximately 25.0942% of total execution time)

with no bining, constant dens 32
Successful octree uses:	120184
Octree fallback uses:	72527
Again, construction time:	59.2701
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 26.2013 seconds.
 Universe::find_cell took 7.97431 seconds in total. (Approximately 30.4348% of total execution time)

after fixing binning
Successful octree uses:	120195
Octree fallback uses:	72620
Again, construction time:	76.3437
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 26.1909 seconds.
 Universe::find_cell took 7.96369 seconds in total. (Approximately 30.4063% of total execution time)

affter seperation (binning)
Successful octree uses:	120149
Octree fallback uses:	72647
Again, construction time:	77.1117
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 26.1907 seconds.
 Universe::find_cell took 7.91638 seconds in total. (Approximately 30.2259% of total execution time)


after bin search
Successful octree uses:	120131
Octree fallback uses:	72324
Again, construction time:	58.6078
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 25.8686 seconds.
 Universe::find_cell took 7.99754 seconds in total. (Approximately 30.9161% of total execution time)

4-12-16
Successful octree uses:	120174
Octree fallback uses:	72944
Again, construction time:	59.5206
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 27.0091 seconds.
 Universe::find_cell took 8.09066 seconds in total. (Approximately 29.9554% of total execution time)

2-2-4-8-16
Successful octree uses:	120140
Octree fallback uses:	72310
Again, construction time:	59.8562
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 26.2316 seconds.
 Universe::find_cell took 7.92618 seconds in total. (Approximately 30.2161% of total execution time)

4-4-4-4-8-8
Successful octree uses:	120151
Octree fallback uses:	72443
Again, construction time:	59.5919
 ============================>    PERFORMANCE STATISTICS    <============================
 Note: performance metrics only count the time spent in the transport loop and currently only support history-based transport.
 The time below is given as total wall clock seconds across all threads (e.g. 1 second of 8 threads executing is counted as 8 seconds)
 Total OpenMC execution time was 26.1838 seconds.
 Universe::find_cell took 7.90228 seconds in total. (Approximately 30.18% of total execution time)
            bool unique = true;
            for(int x : cells) {
                if(x == p.cell) {
                    unique = false;
                    break;
                }
            }

            if(unique) {
                cells.push_back(p.cell);
            }

Algorithm:

construct_octree:
    build root node
    if number of cells in root node is less than target:
        return

    push root node on stack

    while stack not empty:
        pop node from stack and store it in variable "current"

        subdivide current
        for each child
            count cells in child
            if number of cells in child exceeds target
                push child to stack

search_octree:
    current = root
    while true
        if current is leaf
            break
        else
            for each child
                if child contains point
                current = child
                break



proper mt construction algo? multiple threads exec this loop:
    finished = false
    push root to unsplit node stack
    while true
        update info to console if needed

        if finished
            break
        else if unsplit nodes remain
            split nodes
            push split results to unsearched node stack
            wake threads
        else if  unsearched nodes remain
            search node
            if node should be split
                push to unsplit node stack
            else if all other threads are sleeping
                finished = true
            wake threads
        else
            sleep
    end while 

















    older psuedocode:

        if more nodes to search points for
            search points
            if should split
                push to split stack
            else if all other threads are asleep
                mark should exit flag
                wake threads
        else if more unsplit nodes
            split nodes and push to searching stack
            wake threads
        else
            sleep
*/