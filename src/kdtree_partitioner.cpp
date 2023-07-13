#include "openmc/kdtree_partitioner.h"
#include "openmc/error.h"
#include <stack>
#include <queue>
#include <float.h>

namespace openmc {

    const uint32_t KD_TREE_LEAF_FLAG = (0b11 << 30);

    KdTreeUncompressedNode::KdTreeUncompressedNode() : children(nullptr) {}

    bool KdTreeUncompressedNode::is_leaf() const {
        return (children == nullptr);
    }

    bool KdTreeNode::is_leaf() const {
        return ((data & KD_TREE_LEAF_FLAG) == KD_TREE_LEAF_FLAG);
    }

    uint32_t KdTreeNode::index() const {
        uint32_t index = data & ~KD_TREE_LEAF_FLAG;
        return index;
    }


    KdTreeConstructionTask::KdTreeConstructionTask() : node(nullptr), depth(0) {}

    KdTreeConstructionTask::KdTreeConstructionTask(KdTreeUncompressedNode* node, const std::vector<CellPointUncompressed>& points, uint32_t depth) 
                            : node(node), points(points), depth(depth) {}

    uint32_t get_num_unique_cells(std::vector<CellPointUncompressed>& points) {
        std::sort(points.begin(), points.end());

        int prev_cell = -1;
        uint32_t num_unique = 0;
        for(const auto& p : points) {
            if(p.cell != prev_cell) {
                num_unique++;
                prev_cell = p.cell;
            }
        }

        return num_unique;
    }

    void update_best_split_median(const KdTreeConstructionTask& parent_task, KdTreeUncompressedNode best_children[2], std::vector<CellPointUncompressed> best_points[2], float& best_child_vh, int axis) {
        KdTreeUncompressedNode cur_children[2];
        std::vector<CellPointUncompressed> cur_points[2];
        
        // split (code copied from octree)
        float midpoint;
        {

            int i = axis;
            int j = ((i + 1) % 3);
            int k = ((i + 2) % 3);

            midpoint = (parent_task.node->box.max[i] + parent_task.node->box.min[i]) * 0.5;
            const auto& box = parent_task.node->box;

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

            cur_children[0].box = splitted_boxes[1];
            cur_children[1].box = splitted_boxes[0];
        }

        for(const auto& p : parent_task.points) {
            if(cur_children[0].box.contains(p.pos)) {
                cur_points[0].push_back(p);
            } else {
                cur_points[1].push_back(p);
            }
        }

        float cur_vh = 0.0;
        for(int i = 0; i < 2; i++) {
            cur_vh += cur_children[i].box.volume() * get_num_unique_cells(cur_points[i]);
        }

        if(cur_vh < best_child_vh) {
            best_child_vh = cur_vh;

            for(int i = 0; i < 2; i++) {
                best_children[i] = std::move(cur_children[i]);
                best_points[i] = std::move(cur_points[i]);
            }

            parent_task.node->split_axis = axis;
            parent_task.node->split_location = midpoint;
        }
    }

    void make_elems_unique(std::vector<int>& v) {
        #if 0
        std::sort(v.begin(), v.end());
    
        int i = 1, j = 0;
        while(i < v.size()) {
            if(v[i] != v[j]) {
                j++;
                v[j] = v[i];
            }
            i++;
        } 

        v.resize(j);
        #else
        std::set<int> unique_elems;

        for(int x : v) {
            unique_elems.insert(x);
        }

        v.clear();

        for(int x : unique_elems) {
            v.push_back(x);
        }

        #endif
    }

    struct KdTreeBin {
        std::vector<int> unique_cells;
    };

    void consume_bin(const KdTreeBin& bin, std::vector<int>& cells) {
        std::copy(bin.unique_cells.begin(), bin.unique_cells.end(), std::back_insert_iterator(cells));
    }

    int best_picked_split = 0;
    void update_best_split_binning(const KdTreeConstructionTask& parent_task, KdTreeUncompressedNode best_children[2], std::vector<CellPointUncompressed> best_points[2], float& best_child_vh, int axis) {
        const int NUM_KD_TREE_BINS = 64;
        float bin_dim = (parent_task.node->box.max[axis] - parent_task.node->box.min[axis]) / NUM_KD_TREE_BINS;

        KdTreeBin bins[NUM_KD_TREE_BINS];

        for(const auto& p : parent_task.points) {
            int idx = ((p.pos[axis] - parent_task.node->box.min[axis]) / bin_dim);

            if(idx < 0) idx = 0;
            else if(idx >= NUM_KD_TREE_BINS) idx = NUM_KD_TREE_BINS - 1;

            bins[idx].unique_cells.push_back(p.cell);
        }

        for(int i = 0; i < NUM_KD_TREE_BINS; i++) {
            make_elems_unique(bins[i].unique_cells);
        }

        for(int i = 1; i < NUM_KD_TREE_BINS - 1; i++) {
            float split = i * bin_dim + parent_task.node->box.min[axis];
            std::vector<int> cur_cells_list[2]; 

            for(int j = 0; j < i; j++) {
                consume_bin(bins[j], cur_cells_list[0]);
            }

            for(int j = i; j < NUM_KD_TREE_BINS; j++) {
                consume_bin(bins[j], cur_cells_list[1]);
            }

            for(int i = 0; i < 2; i++) {
                make_elems_unique(cur_cells_list[i]);
            }

            // construct AABB
            AABB boxes[2];
            float child_vh = 0.0;
            for(int i = 0; i < 2; i++) {
                boxes[i] = parent_task.node->box;

                auto update = (i == 0 ? &boxes[i].max[axis] : &boxes[i].min[axis]);
                *update = split;

                child_vh += boxes[i].volume() * cur_cells_list[i].size();
            }

            
            if(child_vh < best_child_vh) {
                best_child_vh = child_vh;

                for(int i = 0; i < 2; i++) {
                    best_children[i].box = boxes[i];
                    best_points[i].clear();
                }

                for(const auto& p : parent_task.points) {
                    int idx = (p.pos[axis] < split ? 0 : 1);
                    best_points[idx].push_back(p);
                }

                parent_task.node->split_axis = axis;
                parent_task.node->split_location = split;

                best_picked_split = i;
            }
            
        }

    }


    void make_leaf(KdTreeUncompressedNode& node, std::vector<CellPointUncompressed>& points) {
        node.children = nullptr;

        std::sort(points.begin(), points.end());

        int prev_cell = -1;
        for(const auto& p : points) {
            if(p.cell != prev_cell) {
                node.cells.push_back(p.cell);
                prev_cell = p.cell;
            }
        }
    }

    KdTreePartitioner::KdTreePartitioner(const Universe& univ) : fallback(univ) {
        write_message("Building kd-tree partitioner...", 5);

        Timer construction_timer;
        construction_timer.start();

        const uint32_t MAX_DEPTH = 16;

        const float half_side_length = 130.0;
        bounds.min = vec3(-half_side_length, -half_side_length, -half_side_length);
        bounds.max = vec3( half_side_length,  half_side_length,  half_side_length);

        auto points_in_bounds = binned_point_search<CellPointUncompressed>(univ, fallback, bounds);

        uint32_t num_nodes = 1;
        uint32_t num_leaves = 0;

        NodeAllocator<KdTreeUncompressedNode, 2> node_alloc;

        KdTreeUncompressedNode root;
        root.box = bounds;

        std::stack<KdTreeConstructionTask> unproc_nodes;
        unproc_nodes.emplace(&root, points_in_bounds, 0);
        while(!unproc_nodes.empty()) {
            auto task = std::move(unproc_nodes.top());
            unproc_nodes.pop();

            auto cur = task.node;

            float parent_vh = cur->box.volume() * get_num_unique_cells(task.points);

            KdTreeUncompressedNode best_children[2];
            std::vector<CellPointUncompressed> best_points[2];
            float best_child_vh = FLT_MAX;
            // find best split
            for(int i = 0; i < 3; i++) {
                update_best_split_median(task, best_children, best_points, best_child_vh, i);
            }
            
            // update if needed
            if(0.9 * best_child_vh < parent_vh) {
                num_nodes += 2;

                task.node->children = node_alloc.allocate();
                uint32_t next_depth = task.depth + 1;

                // why am I using a loop for this? I have no idea
                for(int i = 0; i < 2; i++) {
                    task.node->children[i] = std::move(best_children[i]);

                    // make sure depth is less!
                    if(next_depth < MAX_DEPTH) {
                        unproc_nodes.emplace(&task.node->children[i], std::move(best_points[i]), next_depth);
                    } else {
                        make_leaf(task.node->children[i], best_points[i]);
                        num_leaves++;
                    }
                }
            } else {
                make_leaf(*cur, task.points);
                num_leaves++;
            }
        }

        // compress
        nodes.reserve(num_nodes);
        cell_data.reserve(num_leaves);

        std::queue<KdTreeUncompressedNode*> uncomp_nodes;
        uncomp_nodes.push(&root);

        while(!uncomp_nodes.empty()) {
            auto& cur = *uncomp_nodes.front();
            uncomp_nodes.pop();

            // compress
            KdTreeNode comp_node;
            if(cur.is_leaf()) {
                comp_node.data = cell_data.size();
                comp_node.data |= KD_TREE_LEAF_FLAG;

                cell_data.push_back(std::move(cur.cells));
            }  else {
                comp_node.data = (cur.split_axis << 30);
                comp_node.data += nodes.size() + 1 + uncomp_nodes.size();
                comp_node.split = cur.split_location;

                for(int i = 0; i < 2; i++) {
                    uncomp_nodes.push(&cur.children[i]);
                }
            }

            nodes.push_back(comp_node);
        }
        write_message("Kd-tree construction completed in " + std::to_string(construction_timer.elapsed()) + " seconds.");
    }

    KdTreePartitioner::~KdTreePartitioner() {}

    //! Return the list of cells that could contain the given coordinates.
    const std::vector<int32_t>& KdTreePartitioner::get_cells(Position r, Direction u) const {
        if(!bounds.contains(r)) {
            return get_cells_fallback(r, u);
        }

        const KdTreeNode* cur = &nodes[0];
        while(true) {
            uint32_t axis = (cur->data >> 30);
            if(axis == 3) {
                break;
            } else {
                cur = &nodes[cur->index() + uint32_t(r[axis] > cur->split)];
            }
        }

        const auto& cells = cell_data[cur->index()];
        if(cells.empty()) {
            return get_cells_fallback(r, u);
        }

        return cells;
    }

    const std::vector<int32_t>& KdTreePartitioner::get_cells_fallback(Position r, Direction u) const {
        return fallback.get_cells(r, u);
    }


};