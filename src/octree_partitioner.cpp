#include "openmc/octree_partitioner.h"
#include "openmc/random_dist.h"
#include "openmc/timer.h"
#include "openmc/error.h"

#include <omp.h>

#include <algorithm>
#include <stack>
#include <queue>
#include <fstream>

#include <stdlib.h>

namespace openmc {
    const int32_t OMCP_CURRENT_VERSION[3] = {1, 0, 0};

    const float   REFINEMENT_SEARCH_DENSITY          = 0.125;
    const float   REFINEMENT_TIMEOUT                 = 10.0;
    const int32_t REFINEMENT_GRID_RES                = 128; // ideally should be a power of 2

    const int32_t INFORMATION_REFILLING_START_DEPTH_OFFSET = 1;

    enum OMCPStructureType {
        Octree = 1,
        KdTree = 2,               // not yet implemented in OMCP
        ZPlanePartitioner = 3,    // not yet implemented in OMCP
    };

    struct OctreeConstructionTask {
        OctreeUncompressedNode* node;
        std::vector<CellPoint> points;

        OctreeConstructionTask() = default;
        OctreeConstructionTask(OctreeUncompressedNode* n, const std::vector<CellPoint>& p) : node(n), points(p) {}
    };

    OctreeNode::OctreeNode() : data(0) {}

    const uint32_t OCTREE_LEAF_FLAG = (1 << 31);
    void OctreeNode::mark_as_leaf() {
        data |= OCTREE_LEAF_FLAG;
    }

    bool OctreeNode::is_leaf() const {
        return (data & OCTREE_LEAF_FLAG);
    }

    void OctreeNode::store_data(uint32_t data) {
        // remove everything except the flag
        this->data = (this->data & OCTREE_LEAF_FLAG) | data;
    }

    uint32_t OctreeNode::read_data() const {
        return (data & ~OCTREE_LEAF_FLAG);
    }

    OctreeUncompressedNode::OctreeUncompressedNode() : children(nullptr), parent(nullptr), depth(0) {}

    bool OctreeUncompressedNode::is_leaf() const {
        return (children == nullptr);
    }

    void OctreeUncompressedNode::subdivide() {
        AABB resultant_boxes[8];
        resultant_boxes[0] = box;
        for(int i = 0; i < 3; i++) {
            AABB temp_box_buffer[8];

            int next_index = 0;
            for(int idx = 0; idx < (1 << i); idx++) {
                // split on i-th axis
                float midpoint = box.get_center()[i];

                int j = ((i + 1) % 3);
                int k = ((i + 2) % 3);

                AABB splitted_boxes[2] {resultant_boxes[idx], resultant_boxes[idx]};

                splitted_boxes[0].max[i] = midpoint;
                splitted_boxes[1].min[i] = midpoint;

                temp_box_buffer[next_index++] = splitted_boxes[1];
                temp_box_buffer[next_index++] = splitted_boxes[0];
            }

            // move the results to the next splitting stage
            std::copy(temp_box_buffer, temp_box_buffer + (2 << i), resultant_boxes);
        }

        for(int i = 0; i < 8; i++) {
            children[i].box = resultant_boxes[i];
            children[i].depth = depth + 1;
            children[i].parent = this;
        }
    }
    
    bool OctreeUncompressedNode::contains(int cell) const {
        return std::binary_search(cells.begin(), cells.end(), cell);
    }

    void pick_untested_cells(const std::vector<int>& tested_cells, const std::vector<int>& possible_cells, std::vector<int>& untested_cells) {
        untested_cells.clear();
        int next_idx = 0;
        for(int cell : possible_cells) {
            if(next_idx < tested_cells.size() && tested_cells[next_idx] == cell) {
                next_idx++;
            } else {
                untested_cells.push_back(cell);
            }
        }
    }

    void refine_octree_random(
        const Universe& univ,
        const UniversePartitioner& fallback,
        const AABB& bounds,
        const std::vector<OctreeUncompressedNode*>& leaves
    ) {
        const int32_t num_threads = omp_get_max_threads();

        // generate the seeds
        std::vector<uint64_t[2]> rng_node_selec(num_threads);
        std::vector<uint64_t[3]> rng_pos(num_threads);
        for(int i = 0; i < num_threads; i++) {
            rng_node_selec[i][0] = i;
            rng_node_selec[i][1] = i + num_threads;
            for(int j = 0; j < 3; j++) {
                rng_pos[i][j] = i + num_threads * (j + 2);
            }
        }

        using ProbBinT = ProbabilityBin<OctreeUncompressedNode>;
        std::vector<ProbBinT> prob_bin_grid(REFINEMENT_GRID_RES * REFINEMENT_GRID_RES * REFINEMENT_GRID_RES);

        vec3 prob_bin_dim;
        for(int i = 0; i < 3; i++) {
            prob_bin_dim[i] = (bounds.max[i] - bounds.min[i]) / REFINEMENT_GRID_RES;
        }

        for(auto leaf : leaves) {
            int idx = 0;
            for(int i = 0; i < 3; i++) {
                idx = REFINEMENT_GRID_RES * idx + int( (leaf->box.get_center()[i] - bounds.min[i]) / prob_bin_dim[i] );
            }

            omp_lock_t dummy;
            prob_bin_grid[idx].contained_nodes.emplace_back(leaf, dummy);
            for(int cell : leaf->cells) {
                prob_bin_grid[idx].found_cells.push_back(cell);
            }
        }

        for(auto& prob_bin : prob_bin_grid) {
            prob_bin.update_common_cells();
            for(auto& p : prob_bin.contained_nodes) {
                omp_init_lock(&p.second);
            }
        }

        std::vector<float> bin_cdf(prob_bin_grid.size());

        Timer timeout_timer;
        timeout_timer.start();
        int iteration = 0;
        while(timeout_timer.elapsed() < REFINEMENT_TIMEOUT) {
            int num_search_points = (int)(bounds.volume() * REFINEMENT_SEARCH_DENSITY);
            int num_points_searched = 0;

            // first, generate cdf
            float total_cdf = 0.0;
            for(int i = 0; i < prob_bin_grid.size(); i++) {
                total_cdf += prob_bin_grid[i].compute_score();
                bin_cdf[i] = total_cdf;
            }

            for(float& cdf : bin_cdf) {
                cdf /= total_cdf;
            }
            
            #pragma omp parallel for
            for(int tid = 0; tid < num_threads; tid++) {
                std::vector<int> untested_cells, temp_buf;
                untested_cells.reserve(256);
                temp_buf.reserve(256);

                // we don't need a lock for loop iteration logic since nothing bad will happen if race conditions occur
                // maybe we do a few points less or more, but that doesn't matter
                while(num_points_searched < num_search_points) {
                    num_points_searched++;

                    ProbBinT* prob_bin;
                    do {
                        float cdf_val = uniform_distribution(0.0, 1.0, &rng_node_selec[tid][0]);
                        auto cdf_iter = std::upper_bound(bin_cdf.begin(), bin_cdf.end(), cdf_val);
                        size_t bin_idx = std::distance(bin_cdf.begin(), cdf_iter);
                        if(bin_idx == bin_cdf.size()) {
                            bin_idx--;
                        }

                        prob_bin = &prob_bin_grid[bin_idx];
                    } while(prob_bin->contained_nodes.size() == 0);
                    prob_bin->num_searched++;

                    size_t idx = prob_bin->pick_node(&rng_node_selec[tid][1]);
                    if(idx >= prob_bin->contained_nodes.size()) {
                        idx =  prob_bin->contained_nodes.size() - 1;
                    }

                    OctreeUncompressedNode* current = prob_bin->contained_nodes[idx].first;
                    omp_lock_t* cell_lock = &prob_bin->contained_nodes[idx].second;

                    CellPointUncompressed point;
                    for(int i = 0; i < 3; i++) {
                        point.pos[i] = uniform_distribution(current->box.min[i], current->box.max[i], &rng_pos[tid][i]);
                    }

                    omp_set_lock(cell_lock);
                    point.cell = univ.find_cell_for_point(current->cells, point.pos);
                    omp_unset_lock(cell_lock);

                    if(point.cell == -1) {
                        prob_bin->num_found++;

                        omp_set_lock(cell_lock);
                        pick_untested_cells(current->cells, prob_bin->common_cells, untested_cells);
                        omp_unset_lock(cell_lock);

                        point.cell = univ.find_cell_for_point(untested_cells, point.pos);
                        if(point.cell == -1) {
                            point.cell = univ.find_cell_for_point(current->parent->cells, point.pos);

                            if(point.cell == -1) {
                                Direction dummy{0, 0, 1};
                                const auto& possible_cells = fallback.get_cells(point.pos, dummy);

                                pick_untested_cells(untested_cells, possible_cells, temp_buf);

                                omp_set_lock(cell_lock);
                                pick_untested_cells(current->cells, temp_buf, untested_cells);
                                omp_unset_lock(cell_lock);

                                point.cell = univ.find_cell_for_point(untested_cells, point.pos);

                                // very rarely, even the fallback misses a cell
                                // we need to do an exhaustive search or the program will segfault
                                if(point.cell == -1) {
                                    point.cell = univ.find_cell_for_point(univ.cells_, point.pos);
                                }
                            }

                            prob_bin->add_cell(point.cell);
                        }


                        omp_set_lock(cell_lock);
                        // insertion sort
                        current->cells.push_back(point.cell);
                        for(int i = current->cells.size() - 1; i > 0; i--) {
                            if(current->cells[i] < current->cells[i - 1]) {
                                std::swap(current->cells[i], current->cells[i - 1]);
                            }
                        }
                        omp_unset_lock(cell_lock);
                    }
                }
            }
            iteration++;
        }
    }

    OctreePartitioner::OctreePartitioner(const Universe& univ, int target_cells_per_node) : fallback(univ) {
        const float half_side_length = 130.0;
        bounds.min = vec3(-half_side_length, -half_side_length, -half_side_length);
        bounds.max = vec3( half_side_length,  half_side_length,  half_side_length);

        if(univ.cells_.size() <= target_cells_per_node) {
            warning("Universe has only "  + std::to_string(univ.cells_.size()) + 
                    " cells, which is below the target cells per node, which is " 
                    + std::to_string(target_cells_per_node) + " cells. Octree will only consist of root.");
            
            OctreeNode node;
            node.store_data(0);
            node.mark_as_leaf();
            cell_data.push_back(univ.cells_);            
            return;
        }

        write_message("Building octree...", 5);

        Timer construction_timer;
        construction_timer.start();

        auto points_in_bounds = binned_point_search<CellPoint>(univ, fallback, bounds);

        OctreeUncompressedNode root;
        root.box = bounds;
        root.depth = 0;
    
        std::sort(points_in_bounds.begin(), points_in_bounds.end());
        int prev_cell = -1;
        for(const auto& p : points_in_bounds) {
            if(prev_cell != p.get_cell()) {
                root.cells.push_back(p.get_cell());
                prev_cell = p.get_cell();
            }
        }
        root.num_unique_cells = root.cells.size();

        num_nodes = 1;
        num_leaves = 0;

        float depth_vh_mult[] = {
            1.0, 1.0, 1.0, 1.5, 2.5, 4.0,
            6.0, 12.5, 19.0, 32.0, 64.0, 128.0, 
            999.0, 9999.0, 99999.0
        };

        NodeAllocator<OctreeUncompressedNode, 8> node_alloc;
        std::vector<OctreeUncompressedNode*> nodes_to_propagate {&root};
        std::vector<OctreeUncompressedNode*> leaves;

        // this section of code still needs to be multithreaded
        // it can become a bottleneck, espcially with large number of points
        std::queue<OctreeConstructionTask> unprocessed_tasks;
        unprocessed_tasks.emplace(&root, points_in_bounds);
        while(!unprocessed_tasks.empty()) {

            auto cur_task = std::move(unprocessed_tasks.front());
            unprocessed_tasks.pop();

            // subdivide
            cur_task.node->children = node_alloc.allocate();
            cur_task.node->subdivide();

            OctreeConstructionTask child_tasks[8];
            for(int i = 0; i < 8; i++) {
                child_tasks[i].node = &cur_task.node->children[i];
                nodes_to_propagate.push_back(child_tasks[i].node);
            }


            // sort points
            for(const auto& point : cur_task.points) {
                child_tasks[point.get_child_index(cur_task.node->depth)].points.push_back(point);
            }

            float parent_vh = cur_task.node->box.volume() * cur_task.node->num_unique_cells;
            float children_vh = 0.0;

            // post processing (make nodes leaves or push on construction stack)
            bool force_subdiv = false;
            for(int i = 0; i < 8; i++) {
                // count the number of unique cells
                int num_unique_cells = 0;
                int prev_cell = -1;
                for(const auto& p : child_tasks[i].points) {
                    if(p.get_cell() != prev_cell) {
                        num_unique_cells++;
                        prev_cell = p.get_cell();
                    }
                }

                child_tasks[i].node->num_unique_cells = num_unique_cells;

                children_vh += child_tasks[i].node->box.volume() * num_unique_cells;
                if(num_unique_cells > target_cells_per_node) {
                    force_subdiv = true;
                }
            }

            if(force_subdiv || depth_vh_mult[cur_task.node->depth] * children_vh < parent_vh) {
                // continue subdivision on this branch
                num_nodes += 8;
                for(int i = 0; i < 8; i++) {
                    child_tasks[i].node->cells.reserve(child_tasks[i].node->num_unique_cells);
                    prev_cell = -1;
                    for(const auto& p : child_tasks[i].points) {
                        if(p.get_cell() != prev_cell) {
                            child_tasks[i].node->cells.push_back(p.get_cell());
                            prev_cell = p.get_cell();
                        }
                    }

                    unprocessed_tasks.push(std::move(child_tasks[i]));
                }
            } else {
                // terminate subdivision on this branch
                cur_task.node->children = nullptr; 
                leaves.push_back(cur_task.node);
                num_leaves++;
            }


            // free memory
            cur_task.points.clear();
            cur_task.points.shrink_to_fit();
        }

        refine_octree_random(univ, fallback, bounds, leaves);

        // now, build the cells list
        struct NodeComp {
            inline bool operator()(const OctreeUncompressedNode* lhs, const OctreeUncompressedNode* rhs) const {
                return (lhs->depth < rhs->depth);
            }
        };
        std::sort(nodes_to_propagate.rbegin(), nodes_to_propagate.rend(), NodeComp());

        for(auto ptr : nodes_to_propagate) {
            auto& cur = *ptr;

            if(!cur.is_leaf()) {
                std::sort(cur.cells.begin(), cur.cells.end());

                int prev_cell = cur.cells[0];

                int next_idx = 1;
                for(int i = 1; i < cur.cells.size(); i++) {
                    if(cur.cells[i] != prev_cell) {
                        cur.cells[next_idx] = cur.cells[i];
                        next_idx++;
                        prev_cell = cur.cells[i];
                    }
                }
                cur.cells.resize(next_idx);
                cur.cells.shrink_to_fit();
            }

            if(cur.parent) {
                for(int cell : cur.cells) {
                    cur.parent->cells.push_back(cell);
                }
            }
        }

        // a possible way to remove issues with missing information in nodes
        std::stack<OctreeUncompressedNode*> uninfo_refilled_nodes;
        uninfo_refilled_nodes.push(&root);
        while(!uninfo_refilled_nodes.empty()) {
            auto cur = uninfo_refilled_nodes.top();
            uninfo_refilled_nodes.pop();

            bool should_collect_leaves = false;
            for(int i = 0; i < 8; i++) {
                if(cur->children[i].is_leaf()) {
                    should_collect_leaves = true;
                    break;
                }
            }
            
            auto collection_start = cur;
            for(int i = 0; i < INFORMATION_REFILLING_START_DEPTH_OFFSET; i++) {
                if(collection_start->parent) {
                    collection_start = collection_start->parent;
                }
            }
            
            const auto& unique_cells = collection_start->cells;

            for(int i = 0; i < 8; i++) {
                if(cur->children[i].is_leaf()) {
                    int orig_size = cur->children[i].cells.size();

                    auto search_start = cur->children[i].cells.begin();
                    auto search_end   = cur->children[i].cells.begin() + orig_size;

                    cur->children[i].num_original_cells = orig_size;

                    for(int cell : unique_cells) {
                        if(!std::binary_search(search_start, search_end, cell)) {
                            cur->children[i].cells.push_back(cell);
                        }
                    }
                } else {
                    uninfo_refilled_nodes.push(&cur->children[i]);
                }
            }
        }

        // now copy everything to array
        nodes.reserve(num_nodes);
        cell_data.reserve(num_leaves);
        orig_size.reserve(num_leaves);
        std::queue<const OctreeUncompressedNode*> unwritten_nodes;
        unwritten_nodes.push(&root);
        while(!unwritten_nodes.empty()) {
            auto cur = unwritten_nodes.front();
            unwritten_nodes.pop();

            // convert and write
            OctreeNode compressed;

            if(cur->is_leaf()) {
                compressed.store_data(cell_data.size());
                cell_data.push_back(std::move(cur->cells));
                compressed.mark_as_leaf();

                orig_size.push_back(cur->num_original_cells);
            } else {
                compressed.store_data(nodes.size() + 1 + unwritten_nodes.size());
                for(int i = 0; i < 8; i++) {
                    unwritten_nodes.push(&cur->children[i]);
                }            
            }

            nodes.push_back(compressed);
        }

        write_message("Octree construction completed in " + std::to_string(construction_timer.elapsed()) + " seconds.", 5);
    }  

    OctreePartitioner::OctreePartitioner(const Universe& univ, const std::string& path) : fallback(univ) {
        write_message("Reading octree from " + path, 5);

        std::fstream octree_file(path, std::ios::in | std::ios::binary);
        #define READ_BINARY(x) octree_file.read((char*)&x, sizeof(x))

        for(int i = 0; i < 3; i++) {
            int version;
            READ_BINARY(version);
            if(version != OMCP_CURRENT_VERSION[i]) {
                fatal_error("OpenMC cannot read an unsupported OMCP file version! Please note that OpenMC currently cannot read older file versions!");
            }
        }

        OMCPStructureType structure_type;
        READ_BINARY(structure_type);

        if(structure_type != OMCPStructureType::Octree) {
            fatal_error("OpenMC currently only supports octrees for the OMCP file format!");
        }

        std::vector<uint16_t> comp_cell_data;

        READ_BINARY(bounds);

        READ_BINARY(num_nodes);

        int num_leaves = 0;
        nodes.resize(num_nodes);
        for(auto& node : nodes) {
            READ_BINARY(node);
            num_leaves++;
        }
        cell_data.reserve(num_leaves);

        int num_cells;
        READ_BINARY(num_cells);
        comp_cell_data.resize(num_cells);
        for(auto& cell : comp_cell_data) {
            READ_BINARY(cell);
        }

        int next_cell_index = 0;
        for(auto& node : nodes) {
            if(node.is_leaf()) {
                std::vector<int> cells;
                cells.resize(node.read_data());

                for(int& cell : cells) {
                    cell = comp_cell_data.at(next_cell_index);
                    next_cell_index++;
                }

                node.store_data(cell_data.size());
                cell_data.push_back(std::move(cells));
            }
        }

        refill_information();
    }

    OctreePartitioner::~OctreePartitioner() {}

    const vector<int32_t>& OctreePartitioner::get_cells(Position r, Direction u) const {
        if(!bounds.contains(r)) {
            // discount this get_cells call from the stats to only focus on points within the octree
            return get_cells_fallback(r, u);
        }

        vec3 node_dim;
        for(int i = 0; i < 3; i++) {
            node_dim[i] = (bounds.max[i] - bounds.min[i]) * 0.5;
        }

        vec3 center = bounds.get_center();
        auto current = nodes[0];
        while(!current.is_leaf()) {
            // halve the node dim
            int idx = 0;
            for(int i = 0; i < 3; i++) {
                node_dim[i] *= 0.5;
                bool less = (r[i] < center[i]);
                center[i] += node_dim[i] * (less ? -1 : 1);
                idx = 2 * idx + int(less);
            }

            current = nodes[current.read_data() + idx];
        }

        const auto& cells = cell_data[current.read_data()];
        if(cells.empty()) {
            return get_cells_fallback(r, u);
        }

        return cells;
    }

    const vector<int32_t>& OctreePartitioner::get_cells_fallback(Position r, Direction u) const {
        return fallback.get_cells(r, u);
    }

    void OctreePartitioner::export_to_file(const std::string& path) const {
        std::vector<uint16_t> comp_cell_data;
        for(int i = 0; i < cell_data.size(); i++) {
            for(int j = 0; j < orig_size[i]; j++) {
                comp_cell_data.push_back(cell_data[i][j]);
            }
        }
 
        std::fstream octree_file(path, std::ios::out | std::ios::binary);

        #define WRITE_BINARY(x) octree_file.write((char*)&x, sizeof(x))

        for(int i = 0; i < 3; i++) {
            WRITE_BINARY(OMCP_CURRENT_VERSION[i]);
        }
        OMCPStructureType structure_type = OMCPStructureType::Octree;
        WRITE_BINARY(structure_type);

        WRITE_BINARY(bounds);

        WRITE_BINARY(num_nodes);
        for(const auto& raw_node : nodes) {
            auto node = raw_node;
            if(node.is_leaf()) {
                node.store_data(orig_size[node.read_data()]);
            }
            WRITE_BINARY(node);
        }    

        int num_cells = comp_cell_data.size();
        WRITE_BINARY(num_cells);
        for(auto cell : comp_cell_data) {
            WRITE_BINARY(cell);
        }

        write_message("Exported octree to " + path, 5);
    }

    // this method works on the compressed octree, not hte uncompressed one
    // since it is pretty slow (and for some reason creates an octree that has the same failure rate but is slower), 
    // only use it if you are reading from a file
    void OctreePartitioner::refill_information() {
        struct OctreeNodeExtraInfo {
            OctreeNodeExtraInfo() : parent(nullptr), eiparent(nullptr), current(nullptr), eichildren(nullptr), children(nullptr), depth(0) {}

            OctreeNode* parent;
            OctreeNode* current;
            OctreeNode* children;

            OctreeNodeExtraInfo* eiparent;
            OctreeNodeExtraInfo* eichildren;

            std::vector<int> cells;

            uint32_t depth = 0;

            bool is_leaf() const {
                return (children != nullptr);
            }
        };

        // make our nodes easier to work with
        std::vector<OctreeNodeExtraInfo> einodes(nodes.size());
        std::vector<OctreeNodeExtraInfo*> nodes_to_propagate(nodes.size());
        for(int i = 0; i < nodes.size(); i++) {
            OctreeNodeExtraInfo* einode = &einodes[i];
            nodes_to_propagate[i] = einode;

            einode->current = &nodes[i];
            
            if(!nodes[i].is_leaf()) {
                int idx = nodes[i].read_data();
                einode->children = &nodes[idx];
                einode->eichildren = &einodes[idx];

                for(int i = 0; i < 8; i++) {
                    einode->eichildren[i].parent = einode->current;
                    einode->eichildren[i].eiparent = einode;
                    einode->eichildren[i].depth = einode->depth + 1;
                }
            } else {
                einode->cells = std::move(cell_data[nodes[i].read_data()]);
            }
        }

        // propagate all cells forward
        struct DepthComp {
            inline bool operator()(const OctreeNodeExtraInfo* lhs, const OctreeNodeExtraInfo* rhs) const {
                return (lhs->depth < rhs->depth);
            }
        };
        std::sort(nodes_to_propagate.rbegin(), nodes_to_propagate.rend(), DepthComp());
        for(auto ptr : nodes_to_propagate) {
            auto& cur = *ptr;
            if(!cur.is_leaf()) {
                std::sort(cur.cells.begin(), cur.cells.end());

                int prev_cell = cur.cells[0];

                int next_idx = 1;
                for(int i = 1; i < cur.cells.size(); i++) {
                    if(cur.cells[i] != prev_cell) {
                        cur.cells[next_idx] = cur.cells[i];
                        next_idx++;
                        prev_cell = cur.cells[i];
                    }
                }
                cur.cells.resize(next_idx);
                cur.cells.shrink_to_fit();
            }

            if(cur.parent) {
                for(int cell : cur.cells) {
                    cur.eiparent->cells.push_back(cell);
                }
            }
        }

        // now propagate all cells downward
        std::stack<OctreeNodeExtraInfo*> uninfo_refilled_nodes;
        uninfo_refilled_nodes.push(nodes_to_propagate.back());
        while(!uninfo_refilled_nodes.empty()) {
            auto cur = uninfo_refilled_nodes.top();
            uninfo_refilled_nodes.pop();

            bool should_collect_leaves = false;
            for(int i = 0; i < 8; i++) {
                if(cur->children[i].is_leaf()) {
                    should_collect_leaves = true;
                    break;
                }
            }
            
            auto collection_start = cur;
            for(int i = 0; i < INFORMATION_REFILLING_START_DEPTH_OFFSET; i++) {
                if(collection_start->parent) {
                    collection_start = collection_start->eiparent;
                }
            }
            

            const auto& unique_cells = collection_start->cells;

            for(int i = 0; i < 8; i++) {
                if(cur->children[i].is_leaf()) {
                    int orig_size = cur->eichildren[i].cells.size();

                    auto search_start = cur->eichildren[i].cells.begin();
                    auto search_end   = cur->eichildren[i].cells.begin() + orig_size;

                    for(int cell : unique_cells) {
                        if(!std::binary_search(search_start, search_end, cell)) {
                            cur->eichildren[i].cells.push_back(cell);
                        }
                    }
                    // now that downpropagation is done for this node, update cell_data
                    cell_data[cur->eichildren[i].current->read_data()] = std::move(cur->eichildren[i].cells);
                } else {
                    uninfo_refilled_nodes.push(&cur->eichildren[i]);
                }
            }
        }
    }

}