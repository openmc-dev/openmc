#include "openmc/box_tree.h"

#include <cmath>
#include <numeric> // for iota

namespace openmc {

namespace {

//! Minimum number of boxes for a tree to be used. Below it, testing every box
//! is about as fast as searching a tree.
constexpr int MIN_TREE_BOXES = 32;

//! Half of the surface area of a box
double half_area(const BoundingBox& b)
{
  Position d = b.max - b.min;
  return d.x * d.y + d.y * d.z + d.z * d.x;
}

//! Builds a tree top-down using the binned surface area heuristic
class Builder {
public:
  static constexpr int N_BINS = 16;
  static constexpr int MAX_LEAF = 4;
  static constexpr int MAX_SAH_DEPTH = 40;

  Builder(const vector<BoundingBox>& boxes) : boxes_(boxes)
  {
    centers_.reserve(boxes.size());
    for (const auto& b : boxes)
      centers_.push_back(0.5 * (b.min + b.max));
  }

  //! Build the subtree over items [begin, end) of order, appending nodes and
  //! moving the items of each leaf to the end of leaf_order
  template<typename Node>
  void build(vector<Node>& nodes, vector<int32_t>& order, int32_t begin,
    int32_t end, int depth, vector<int32_t>& leaf_order)
  {
    int32_t i_node = nodes.size();
    nodes.push_back({});
    BoundingBox box = BoundingBox::inverted();
    BoundingBox center_box = BoundingBox::inverted();
    for (int32_t k = begin; k < end; ++k) {
      box |= boxes_[order[k]];
      center_box.expand(centers_[order[k]]);
    }
    nodes[i_node].box = box;
    int32_t n = end - begin;

    // Find the split along the axis and bin boundary with the lowest cost
    int best_axis = -1;
    int best_split = 0;
    double best_cost = INFTY;
    if (n > 1 && depth < MAX_SAH_DEPTH) {
      for (int axis = 0; axis < 3; ++axis) {
        double lo = center_box.min[axis];
        double extent = center_box.max[axis] - lo;
        if (!(extent > 0.0))
          continue;
        array<BoundingBox, N_BINS> bin_box;
        array<int, N_BINS> bin_count {};
        bin_box.fill(BoundingBox::inverted());
        for (int32_t k = begin; k < end; ++k) {
          int b = bin(centers_[order[k]][axis], lo, extent);
          bin_box[b] |= boxes_[order[k]];
          ++bin_count[b];
        }
        // Sweep from the right to get the cost of the right side of each
        // split, then from the left
        array<double, N_BINS> right_cost;
        BoundingBox acc = BoundingBox::inverted();
        int count = 0;
        for (int b = N_BINS - 1; b > 0; --b) {
          acc |= bin_box[b];
          count += bin_count[b];
          right_cost[b] = count > 0 ? half_area(acc) * count : 0.0;
        }
        acc = BoundingBox::inverted();
        count = 0;
        for (int b = 1; b < N_BINS; ++b) {
          acc |= bin_box[b - 1];
          count += bin_count[b - 1];
          if (count == 0 || count == n)
            continue;
          double cost = half_area(acc) * count + right_cost[b];
          if (cost < best_cost) {
            best_cost = cost;
            best_axis = axis;
            best_split = b;
          }
        }
      }
    }

    // Make a leaf if splitting does not reduce the expected cost of testing
    // the items, relative to visiting both children
    double leaf_cost = half_area(box) * n;
    bool split_sah = best_axis >= 0 && best_cost < leaf_cost;
    if (n <= MAX_LEAF && !split_sah) {
      nodes[i_node].first = leaf_order.size();
      nodes[i_node].count = n;
      for (int32_t k = begin; k < end; ++k)
        leaf_order.push_back(order[k]);
      return;
    }

    int32_t mid;
    if (split_sah) {
      double lo = center_box.min[best_axis];
      double extent = center_box.max[best_axis] - lo;
      auto it = std::partition(
        order.begin() + begin, order.begin() + end, [&](int32_t item) {
          return bin(centers_[item][best_axis], lo, extent) < best_split;
        });
      mid = it - order.begin();
    } else {
      // Split at the median along the longest axis of the centers
      Position d = center_box.max - center_box.min;
      int axis = d.x >= d.y && d.x >= d.z ? 0 : (d.y >= d.z ? 1 : 2);
      mid = begin + n / 2;
      std::nth_element(order.begin() + begin, order.begin() + mid,
        order.begin() + end, [&](int32_t a, int32_t b) {
          return centers_[a][axis] < centers_[b][axis];
        });
    }

    nodes[i_node].count = 0;
    build(nodes, order, begin, mid, depth + 1, leaf_order);
    nodes[i_node].first = nodes.size();
    build(nodes, order, mid, end, depth + 1, leaf_order);
  }

private:
  static int bin(double x, double lo, double extent)
  {
    int b = static_cast<int>(N_BINS * (x - lo) / extent);
    return std::min(std::max(b, 0), N_BINS - 1);
  }

  const vector<BoundingBox>& boxes_;
  vector<Position> centers_;
};

} // namespace

bool use_box_tree(const vector<BoundingBox>& boxes)
{
  if (boxes.size() < MIN_TREE_BOXES)
    return false;

  // A tree only helps if each point is in few of the boxes. The total volume
  // of the boxes divided by the volume of the box around all of them is the
  // average number of boxes containing a point in it, which must be at most
  // an eighth of the number of boxes.
  auto volume = [](const BoundingBox& b) {
    return (b.max.x - b.min.x) * (b.max.y - b.min.y) * (b.max.z - b.min.z);
  };
  BoundingBox all = BoundingBox::inverted();
  double total = 0.0;
  for (const auto& b : boxes) {
    all |= b;
    total += volume(b);
  }
  return 8.0 * total <= boxes.size() * volume(all);
}

BoxTree::BoxTree(const vector<BoundingBox>& boxes)
{
  if (boxes.empty())
    return;

  // Enlarge the boxes to tolerate roundoff in positions on their surfaces
  vector<BoundingBox> padded;
  padded.reserve(boxes.size());
  for (BoundingBox b : boxes) {
    for (int i = 0; i < 3; ++i) {
      double pad =
        1e-6 + 1e-9 * std::max(std::abs(b.min[i]), std::abs(b.max[i]));
      b.min[i] -= pad;
      b.max[i] += pad;
    }
    padded.push_back(b);
  }

  vector<int32_t> order(boxes.size());
  std::iota(order.begin(), order.end(), 0);
  Builder builder(padded);
  builder.build(nodes_, order, 0, order.size(), 0, items_);
  for (int32_t item : items_)
    item_boxes_.push_back(padded[item]);
}

} // namespace openmc
