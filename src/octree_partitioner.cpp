#include "openmc/octree_partitioner.h"
#include "openmc/aabb.h"
#include "openmc/random_dist.h"
#include "openmc/timer.h"

#include <stack>
#include <thread>
#include <assert.h>
#include <mutex>

//#define CXX_THREAD_MT_CONSTRUCTION
#ifdef CXX_THREAD_MT_CONSTRUCTION
#include <condition_variable>
#endif

#define OCTREE_DEBUG_INFO

namespace openmc {


int octree_use_count = 0;
int fallback_use_count = 0;

#ifdef OCTREE_DEBUG_INFO
std::mutex octree_console_mutex;
#endif

OctreeNode::OctreeNode() : children(nullptr), depth(0) {}

bool OctreeNode::is_leaf() const {
    return (children == nullptr);
}
 
void find_cells_in_box(const Universe& univ, const OctreeNode& parent, OctreeNode& child) {
    const double NUM_CELL_SEARCH_POINT_DENSITY = 1024.0;
    const int num_search_points = (int)(child.box.volume() * NUM_CELL_SEARCH_POINT_DENSITY);

    // get some numbers that are unique across all invocations
    uint64_t octree_prng_seed[3] = {(uint64_t)&child, (uint64_t)child.id, (uint64_t)((child.id * 3244230925 + 432534) % 436786)};

    #ifdef OCTREE_DEBUG_INFO
    const int POINTS_LOGGING_THRESHOLD = 16384; 
    const int POINTS_LOGGING_INCREMENT = 5000;
    if(num_search_points > POINTS_LOGGING_THRESHOLD) {
        octree_console_mutex.lock();
        std::cout << "Note: searching node " << child.id << " for " << num_search_points << " points.\n";
        octree_console_mutex.unlock();
    }

    Timer t;
    t.start();
    #endif

    std::vector<bool> skip_cell(parent.cells.size());
    fill(skip_cell.begin(), skip_cell.end(), false); // maybe it initializes to false? but just making sure

    //std::vector<int> temp_buf;
    //temp_buf.reserve(parent.cells.size());

    child.cells.reserve(parent.cells.size() * 2);
    for(int i = 0; i < num_search_points; i++) {
        #ifdef OCTREE_DEBUG_INFO
        if(num_search_points > POINTS_LOGGING_THRESHOLD && i % POINTS_LOGGING_INCREMENT == 0 && i != 0) {
            octree_console_mutex.lock();
            std::cout << "Node " << child.id << ":\tcurrently searched " << i << " out of " << num_search_points << " points (" << (100.0 * i) / num_search_points 
                      << "%) and found " << child.cells.size() << " unique cells. ";
            std::cout << "ETA is " << ((float)num_search_points / i * t.elapsed()) - t.elapsed() << " more seconds for this node.\n";
            std::cout.flush(); // you may want to comment this out
            octree_console_mutex.unlock();
        }
        #endif

        vec3 rand_pos;

        // gen random point
        // #pragma omp critical
        {
            rand_pos.x = uniform_distribution(child.box.min.x, child.box.max.x, &octree_prng_seed[0]);
            rand_pos.y = uniform_distribution(child.box.min.y, child.box.max.y, &octree_prng_seed[1]);
            rand_pos.z = uniform_distribution(child.box.min.z, child.box.max.z, &octree_prng_seed[2]);
        }

        univ.find_cell_in_list(parent.cells, child.cells, skip_cell, rand_pos);

        // #pragma omp critical
        /*
        {
            for(int cell : temp_buf) {
                child.cells.push_back(cell);
            }
        }
        */
    }

    

#if 0
    if(child.cells.empty()) {
        return;
    }

    // make sure that the vector only contains unique elements
    std::sort(child.cells.begin(), child.cells.end());
    int i = 1, j = 1;
    int prev_value = child.cells[0];
    for(; j < child.cells.size(); j++) {
        if(child.cells[j] == prev_value) {
            continue;
        } else {
            child.cells[i] = child.cells[j];
            prev_value = child.cells[j];
            i++;
        }
    }
    child.cells.resize(i);
#else
    std::set<int> unique_cells;
    for(int cell : child.cells) {
        unique_cells.emplace(cell);
    }
    child.cells.clear();
    for(int cell : unique_cells) {
        child.cells.push_back(cell);
    }
#endif
}

void construction_reporter_thread(int max_tree_size, int* num_nodes, bool* currently_constructing) {
    using namespace std::chrono_literals;
    std::this_thread::sleep_for(500ms); // wait for construction to get going

    int last_num_nodes = *num_nodes;
    while(*currently_constructing) {
        if(*num_nodes != last_num_nodes) {
            std::cout << "Tree size is " << *num_nodes << "\t(" << (100.0 * *num_nodes) / max_tree_size << "% of max tree size)\n";
            last_num_nodes = *num_nodes;
        }
        std::this_thread::sleep_for(1000ms);
    }
}

#ifdef CXX_THREAD_MT_CONSTRUCTION 
void process_child(
    const Universe* univ, const OctreeNode* parent, OctreeNode
  int successful_use_count;
  int fallback_use_count;ndition_variable* launcher_waker
) {
    alive_threads++;
    find_cells_in_box(*univ, *parent, *child);
    if(child->cells.size() > target_cells_per_node) {
        push_mutex->lock();
        unproc_nodes->push(child);
        push_mutex->unlock();
    }
    alive_threads--;
    launcher_waker->notify_all();
}
#endif

Timer t;
OctreePartitioner::OctreePartitioner(const Universe& univ, int target_cells_per_node, int max_depth) : fallback(univ) {
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

    root.box = rootbox;
    root.cells = univ.cells_;
    root.depth = 0;
    root.id = 0;

    std::sort(root.cells.begin(), root.cells.end());

    if(root.cells.size() <= target_cells_per_node) {
        std::cout << "Root has only " << root.cells.size() << " cells, which is below the target cells per node, " << target_cells_per_node << ". Octree will have depth of only 1.\n";
        return;
    }

    #ifdef CXX_THREAD_MT_CONSTRUCTION
    int alive_threads = 0;
    int available_threads = 15; // TODO: fetch this from the OS
    std::mutex launch_mutex;
    std::condition_variable launcher_waker;
    std::vector<std::thread*> searcher_threads;
    #endif

    int next_id = 1;
    int num_nodes = 1;
    bool currently_constructing = true;
    std::thread reporter_thread(construction_reporter_thread, max_tree_size, &num_nodes, &currently_constructing);

    std::mutex push_mutex;

    std::stack<OctreeNode*> unproc_nodes;
    unproc_nodes.push(&root);
    while(!unproc_nodes.empty()) {
        //std::cout << "Processing node...\n";
        auto& current = *unproc_nodes.top();
        unproc_nodes.pop();

        if(current.depth == max_depth) {
            continue;
        }

        // subdivide
        current.children = new OctreeNode[8]; // we'll probably want to use more efficient allocation for this
        num_nodes += 8;
        
        // to write clean code for the subidivision process, I don't hardcode everything out and instead use loops that repeadetly split the box on an axis
        std::vector<AABB> resultant_boxes;
        resultant_boxes.push_back(current.box);
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
            current.children[i].box = resultant_boxes[i];
            current.children[i].depth = current.depth + 1;
            current.children[i].id = next_id;
            next_id++;
        }

        // for each child
        #ifdef CXX_THREAD_MT_CONSTRUCTION
        #error C++ Thread MT is not finished yet
        searcher_threads.clear();
        for(int i = 0; i < 8; i++) {
            bool had_to_wait = false;
            
            while(alive_threads < available_threads) {
                launcher_waker.wait();
            }

            launch_mutex.lock();
            searcher_threads.push_back(new std::thread(process_child, 
                &univ, &current, &current.children[i], &unproc_nodes,
                &push_mutex, target_cells_per_node, &alive_threads, &launcher_waker
            ));
            launch_mutex.unlock();
        }

        for(auto worker : searcher_threads) {
            worker->join();
        }
        #else
        #pragma omp parallel for
        for(int i = 0; i < 8; i++) {
            find_cells_in_box(univ, current, current.children[i]);
            if(current.children[i].cells.size() > target_cells_per_node) {
                push_mutex.lock();
                unproc_nodes.push(&current.children[i]);
                push_mutex.unlock();
            }
        }
        #endif

        // free up memory
        current.cells.clear();
    }    

    t.stop();
    std::cout << "Construction took " << t.elapsed() << " seconds.\n";

    currently_constructing = false;
    reporter_thread.join();

}

OctreePartitioner::~OctreePartitioner() {
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


    if(!root.box.contains(r)) {
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
        bool inside = false;
        for(int i = 0; i < 8; i++) {
            const OctreeNode& child = current->children[i];
            if(child.box.contains(r)) {
                current = &child;
                inside = true;
                break;
            }
        }
        
        if(!inside) {
            #ifdef OCTREE_DEBUG_INFO
            octree_console_mutex.lock();
            #endif
            std::cout << "FATAL ERROR: particle at " << r << " is not contained by any child! It is likely that floating point precision is too low!" << std::endl;
            #ifdef OCTREE_DEBUG_INFO
            octree_console_mutex.unlock();
            #endif
            abort();
        }
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

*/