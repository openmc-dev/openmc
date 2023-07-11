#include "openmc/kdtree_partitioner.h"
#include <stack>
#include <float.h>

namespace openmc {

    KdTreeNode::KdTreeNode() : children(nullptr) {}

    bool KdTreeNode::is_leaf() const {
        return (children == nullptr);
    }

    KdTreeConstructionTask::KdTreeConstructionTask() : node(nullptr), depth(0) {}

    KdTreeConstructionTask::KdTreeConstructionTask(KdTreeNode* node, const std::vector<CellPointUncompressed>& points, uint32_t depth) 
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

    void update_best_split(const KdTreeConstructionTask& parent_task, KdTreeNode best_children[2], std::vector<CellPointUncompressed> best_points[2], float& best_child_vh, int axis) {
        KdTreeNode cur_children[2];
        std::vector<CellPointUncompressed> cur_points[2];
        
        // split (code copied from octree)
        {

            int i = axis;
            int j = ((i + 1) % 3);
            int k = ((i + 2) % 3);

            float midpoint = (parent_task.node->box.max[i] + parent_task.node->box.min[i]) * 0.5;
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

            cur_children[0].box = splitted_boxes[0];
            cur_children[1].box = splitted_boxes[1];
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
        }
    }

    void make_leaf(KdTreeNode& node, std::vector<CellPointUncompressed>& points) {
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
        std::cout << "====================KDTREE CONSTRUCTION===================\n";
        const uint32_t MAX_DEPTH = 16;

        const float half_side_length = 130.0;
        bounds.min = vec3(-half_side_length, -half_side_length, -half_side_length);
        bounds.max = vec3( half_side_length,  half_side_length,  half_side_length);
        root.box = bounds;

        auto points_in_bounds = binned_point_search<CellPointUncompressed>(univ, fallback, bounds);
        std::cout << "Done searching points. Beginning organization of points..." << std::endl;

        std::stack<KdTreeConstructionTask> unproc_nodes;
        unproc_nodes.emplace(&root, points_in_bounds, 0);

        Timer construction_timer;
        construction_timer.start();
        while(!unproc_nodes.empty()) {
            auto task = std::move(unproc_nodes.top());
            unproc_nodes.pop();

            auto cur = task.node;

            float parent_vh = cur->box.volume() * get_num_unique_cells(task.points);
            

            KdTreeNode best_children[2];
            std::vector<CellPointUncompressed> best_points[2];
            float best_child_vh = FLT_MAX;
            // find best split
            for(int i = 0; i < 3; i++) {
                update_best_split(task, best_children, best_points, best_child_vh, i);
            }
            
            // update if needed
            if(best_child_vh < parent_vh) {
                task.node->children = new KdTreeNode[2];
                uint32_t next_depth = task.depth + 1;

                // why am I using a loop for this? I have no idea
                for(int i = 0; i < 2; i++) {
                    task.node->children[i] = std::move(best_children[i]);

                    // make sure depth is less!
                    if(next_depth < MAX_DEPTH) {
                        unproc_nodes.emplace(&task.node->children[i], std::move(best_points[i]), next_depth);
                    } else {
                        make_leaf(task.node->children[i], best_points[i]);
                    }
                }
            } else {
                make_leaf(*cur, task.points);
            }
        }
        construction_timer.stop();
        std::cout << "Kd tree construction took " << construction_timer.elapsed() << " seconds." << std::endl;
    }

    // statistics
    uint32_t kd_tree_use_count = 0;
    uint32_t kd_tree_fallback_use_count = 0;
    uint32_t kd_tree_oob_use_count = 0;
    KdTreePartitioner::~KdTreePartitioner() {
        // memory leak for now

        std::cout << "=========KD TREE USE STATS==========\n";
        std::cout << "Kd tree use count      :\t" << kd_tree_use_count << '\n';
        std::cout << "Fallback use count     :\t" << kd_tree_fallback_use_count << '\n';
        std::cout << "Out of bounds use count:\t" << kd_tree_oob_use_count << '\n';
        std::cout.flush();
    }

    //! Return the list of cells that could contain the given coordinates.
    const std::vector<int32_t>& KdTreePartitioner::get_cells(Position r, Direction u) const {
        if(!bounds.contains(r)) {
            kd_tree_oob_use_count++;
            return get_cells_fallback(r, u);
        }
        kd_tree_use_count++;

        const KdTreeNode* cur = &root;
        while(!cur->is_leaf()) {
            if(cur->children[0].box.contains(r)) {
                cur = &cur->children[0];
            } else {
                cur = &cur->children[1];
            }
        }

        const auto& cells = cur->cells;
        if(cells.empty()) {
            return get_cells_fallback(r, u);
        }

        return cells;
    }

    const std::vector<int32_t>& KdTreePartitioner::get_cells_fallback(Position r, Direction u) const {
        kd_tree_fallback_use_count++;
        return fallback.get_cells(r, u);
    }


};