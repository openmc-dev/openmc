#include "openmc/octree_partitioner.h"
#include "openmc/aabb.h"
#include "openmc/random_dist.h"
#include "openmc/timer.h"

#include <stack>
#include <queue>
#include <thread>
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

#ifdef OCTREE_PARTITIONER_BLOCK_ALLOCATION
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
#endif

Timer t;
const double NUM_CELL_SEARCH_POINT_DENSITY = 32.0;
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

    const int num_search_points = (int)(bounds.volume() * NUM_CELL_SEARCH_POINT_DENSITY);
    std::vector<CellPoint> points_in_bounds(num_search_points);

    omp_lock_t progess_display_lock;
    omp_init_lock(&progess_display_lock);

    int total_points_searched = 0;
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


    std::cout << "Done searching for points! Beginning organization of points in octree...\n";

    struct OctreeConstructionTask {
        OctreeNode* node;
        AABB box;
        std::vector<CellPoint> points;

        OctreeConstructionTask() = default;
        OctreeConstructionTask(OctreeNode* n, const AABB& b, const std::vector<CellPoint>& p) : node(n), box(b), points(p) {}
    };

    std::stack<OctreeConstructionTask> oct_unproc;
    oct_unproc.emplace(&root, bounds, points_in_bounds);
    while(!oct_unproc.empty()) {
        auto cur_task = oct_unproc.top();
        oct_unproc.pop();

        // subdivide
        #ifdef OCTREE_PARTITIONER_BLOCK_ALLOCATION
        cur_task.node->children = node_alloc.allocate();
        #else
        cur_task.node->children = new OctreeNode[8]; 
        #endif

        auto boxes = cur_task.node->subdivide(cur_task.box);

        OctreeConstructionTask child_tasks[8];
        for(int i = 0; i < 8; i++) {
            child_tasks[i].box = boxes[i];
            child_tasks[i].node = &cur_task.node->children[i];
        }

        for(const auto& point : cur_task.points) {
            child_tasks[cur_task.node->get_containing_child_index(point.pos)].points.push_back(point);
        }

        for(int i = 0; i < 8; i++) {
            sstd:;set<int> unique_cells;

            for(const auto& point : child_tasks[i].points) {
                unique_cells.emplace(point.cell);
            }

            if(unique_cells.size() > target_cells_per_node) {
                oct_unproc.push(child_tasks[i]);
            } else {
                child_tasks[i].node->cells.resize(unique_cells.size());
                for(int cell : unique_cells) {
                    child_tasks[i].node->cells.push_back(cell);
                }
            }
        }
    }

    std::cout << "Done sorting points into octree!\n";









#if 0
    int next_id = 1;
    num_nodes = 1;

    omp_lock_t push_lock;

    std::stack<std::pair<OctreeNode*, AABB>> unproc_nodes;
    unproc_nodes.emplace(&root, rootbox);
    while(!unproc_nodes.empty()) {
        //std::cout << "Processing node...\n";
        auto p = unproc_nodes.top();
        auto& current = *p.first;
        auto box = p.second;
        unproc_nodes.pop();

        if(current.depth == max_depth) {
            continue;
        }

        // subdivide
        #ifdef OCTREE_PARTITIONER_BLOCK_ALLOCATION
        current.children = node_alloc.allocate();
        #else
        current.children = new OctreeNode[8]; 
        #endif
        num_nodes += 8;
        
        // to write clean code for the subidivision process, I don't hardcode everything out and instead use loops that repeadetly split the box on an axis
        std::vector<AABB> resultant_boxes;
        resultant_boxes.push_back(box);
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
            current.children[i].center = resultant_boxes[i].get_center();
            current.children[i].depth = current.depth + 1;
            current.children[i].id = next_id;
            next_id++;
        }


        #pragma omp parallel for
        for(int i = 0; i < 8; i++) {
            find_cells_in_box(univ, current, current.children[i], resultant_boxes[i]);
            if(current.children[i].cells.size() > target_cells_per_node) {
                omp_set_lock(&push_lock);
                unproc_nodes.emplace(&current.children[i], resultant_boxes[i]);
                omp_unset_lock(&push_lock);
            }
        }
    }    
#endif

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
    #ifndef OCTREE_PARTITIONER_BLOCK_ALLOCATION
    std::stack<OctreeNode*> unprocessed_nodes;
    unprocessed_nodes.push(&root);

    std::vector<OctreeNode*> unfreed_pointers;
    while(!unprocessed_nodes.empty()) {
        auto& current = *unprocessed_nodes.top();
        unprocessed_nodes.pop();

        if(current.is_leaf()) {
            unfreed_pointers.push_back(current.children);
            for(int i = 0; i < 8; i++) {
                unprocessed_nodes.push(&current.children[i]);
            }
        }
    }

    for(auto ptr : unfreed_pointers) {
        delete[] ptr;
    }


    #endif
    // else, the destructor of the node alloc will take care of things for us

    std::cout << "Successful octree uses:\t" << octree_use_count << '\n';
    std::cout << "Octree fallback uses:\t" << fallback_use_count << '\n';
    std::cout << "Again, construction time:\t" << t.elapsed() << '\n';
}

//#define FALLBACK_OPTIMIZATIONS

const vector<int32_t>& OctreePartitioner::get_cells(Position r, Direction u) const {
    octree_use_count++;

    // completely ignore the direction argument as we wont need that
    
    //#ifdef OCTREE_DEBUG_INFO
    //octree_console_mutex.lock();
    //std::cout << "Searching for particle at " << r << std::endl;
    //octree_console_mutex.unlock();
    //#endif


    if(!bounds.contains(r)) {
        #ifdef FALLBACK_OPTIMIZATIONS
        return get_cells_fallback(r, u);
        #else
        // I don't want this to be counted in the stats
        fallback_use_count--; 
        octree_use_count--;
        // root.cells is empty
        return root.cells;
        #endif
    }


    const OctreeNode* current = &root;

    while(!current->is_leaf()) {
        current = &current->get_containing_child(r);
    }

    #ifdef FALLBACK_OPTIMIZATIONS
    if(current->cells.empty()) {
        return get_cells_fallback(r, u);
    }
    #endif

    return current->cells;
}

const vector<int32_t>& OctreePartitioner::get_cells_fallback(Position r, Direction u) const {
    fallback_use_count++;
    //#ifdef OCTREE_DEBUG_INFO
    //octree_console_mutex.lock();
    //std::cout << "Falling back to z-plane partitioner for particle at " << r << std::endl;
    //octree_console_mutex.unlock();
    //#endif
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