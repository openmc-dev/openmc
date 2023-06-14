#include "openmc/octree_partitioner.h"
#include "openmc/aabb.h"
#include "openmc/random_dist.h"

#include <stack>
#include <assert.h>

namespace openmc {

OctreeNode::OctreeNode() : children(nullptr) {}

bool OctreeNode::is_leaf() const {
    return (children == nullptr);
}

static uint64_t octree_prng_seed[3] = {123, 456, 789};
void find_cells_in_box(const Universe& univ, const OctreeNode& parent, OctreeNode& child) {
    const int NUM_CELL_SEARCH_POINTS = 100;

    for(int i = 0; i < NUM_CELL_SEARCH_POINTS; i++) {
        // gen random point
        vec3 rand_pos;
        rand_pos.x = uniform_distribution(child.box.min.x, child.box.max.x, &octree_prng_seed[0]);
        rand_pos.y = uniform_distribution(child.box.min.y, child.box.max.y, &octree_prng_seed[1]);
        rand_pos.z = uniform_distribution(child.box.min.z, child.box.max.z, &octree_prng_seed[2]);

        univ.find_cell_in_list(parent.cells, child.cells, rand_pos);
    }

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
}

OctreePartitioner::OctreePartitioner(const Universe& univ, int target_cells_per_node) {
    // octrees don't have to have cube nodes/leaves, so we don't have to convert the source AABB into a square
    auto urootbox = univ.bounding_box();
    AABB rootbox(vec3(urootbox.xmin, urootbox.ymin, urootbox.zmin), vec3(urootbox.xmax, urootbox.ymax, urootbox.zmax));

    root.box = rootbox;
    root.cells = univ.cells_;

    if(root.cells.size() <= target_cells_per_node) {
        return;
    }

    std::stack<OctreeNode*> unprocessed_nodes;
    unprocessed_nodes.push(&root);
    while(!unprocessed_nodes.empty()) {
        auto& current = *unprocessed_nodes.top();
        unprocessed_nodes.pop();

        // subdivide
        current.children = new OctreeNode[8]; // we'll probably want to use more efficient allocation for this
        
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
        }

        // for each child
        for(int i = 0; i < 8; i++) {
            find_cells_in_box(univ, current, current.children[i]);
            if(current.children[i].cells.size() > target_cells_per_node) {
                unprocessed_nodes.push(&current.children[i]);
            }
        }
    }    

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
}


const vector<int32_t>& OctreePartitioner::get_cells(Position r, Direction u) const {
    // completely ignore the direction argument as we wont need that
    const OctreeNode* current = &root;

    while(!current->is_leaf()) {
        for(int i = 0; i < 8; i++) {
            const OctreeNode& child = current->children[i];
            if(child.box.contains(r)) {
                current = &child;
                break;
            }
        }
    }

    return current->cells;
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