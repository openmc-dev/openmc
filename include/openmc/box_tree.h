#ifndef OPENMC_BOX_TREE_H
#define OPENMC_BOX_TREE_H

#include <algorithm> // for max, min
#include <cstdint>
#include <limits>

#include "openmc/array.h"
#include "openmc/bounding_box.h"
#include "openmc/constants.h"
#include "openmc/position.h"
#include "openmc/vector.h"

namespace openmc {

//! Whether a box is finite and not empty
inline bool is_bounded(const BoundingBox& b)
{
  for (int i = 0; i < 3; ++i) {
    if (!(b.min[i] > -INFTY && b.max[i] < INFTY && b.min[i] <= b.max[i]))
      return false;
  }
  return true;
}

//==============================================================================
//! Bounding volume hierarchy over a set of axis-aligned boxes.
//!
//! The tree is built top-down by splitting the items at each node with the
//! binned surface area heuristic. It answers which boxes may contain a point
//! and which boxes a ray may enter.
//==============================================================================

class BoxTree {
public:
  BoxTree() = default;

  //! Build a tree over finite boxes. Items are identified by their position
  //! in the vector. Boxes are enlarged slightly so that a point on the surface
  //! of a box, up to roundoff, is considered to be inside of it.
  explicit BoxTree(const vector<BoundingBox>& boxes);

  bool empty() const { return nodes_.empty(); }

  //! Call f(item) for each item whose box contains the point r, stopping as
  //! soon as f returns true.
  //! \return Whether f returned true
  template<typename F>
  bool any_containing(Position r, F&& f) const;

  //! Call f(item, t) for items whose box the ray from r along u enters at a
  //! distance t <= t_max, visiting nearer boxes first. The callable may lower
  //! t_max to prune the remaining search.
  template<typename F>
  void visit_ray(Position r, Direction u, double& t_max, F&& f) const;

private:
  //! A node of the tree. The first child of an interior node immediately
  //! follows it and the index of the second child is stored in first.
  struct Node {
    BoundingBox box;
    int32_t first; //!< First item of a leaf or second child of an interior node
    int32_t count; //!< Number of items of a leaf (0 for an interior node)
  };

  static constexpr int MAX_DEPTH = 64;

  //! Distance along the ray at which it enters the box, or infinity if it
  //! misses. Infinity is larger than any t_max, which is at most INFTY, so a
  //! box that is missed is never visited.
  static double enter(
    const BoundingBox& b, Position r, Position inv_u, double t_max)
  {
    // Comparisons are ordered so that a NaN from a ray parallel to and lying
    // in a plane of the box leaves the bounds unchanged
    double t0 = 0.0;
    double t1 = t_max;
    for (int i = 0; i < 3; ++i) {
      double ta = (b.min[i] - r[i]) * inv_u[i];
      double tb = (b.max[i] - r[i]) * inv_u[i];
      t0 = std::max(t0, std::min(ta, tb));
      t1 = std::min(t1, std::max(ta, tb));
    }
    return t0 <= t1 ? t0 : std::numeric_limits<double>::infinity();
  }

  static bool contains(const BoundingBox& b, Position r)
  {
    return r.x >= b.min.x && r.x <= b.max.x && r.y >= b.min.y &&
           r.y <= b.max.y && r.z >= b.min.z && r.z <= b.max.z;
  }

  vector<Node> nodes_;
  vector<int32_t> items_;          //!< Items in leaf order
  vector<BoundingBox> item_boxes_; //!< Boxes of the items in leaf order
};

//==============================================================================
// Template implementations
//==============================================================================

template<typename F>
bool BoxTree::any_containing(Position r, F&& f) const
{
  if (nodes_.empty())
    return false;
  array<int32_t, MAX_DEPTH> stack;
  int n = 0;
  int32_t i = 0;
  while (true) {
    const Node& node = nodes_[i];
    if (contains(node.box, r)) {
      if (node.count == 0) {
        stack[n++] = node.first;
        i = i + 1;
        continue;
      }
      for (int32_t k = node.first; k < node.first + node.count; ++k) {
        if (contains(item_boxes_[k], r) && f(items_[k]))
          return true;
      }
    }
    if (n == 0)
      return false;
    i = stack[--n];
  }
}

template<typename F>
void BoxTree::visit_ray(Position r, Direction u, double& t_max, F&& f) const
{
  if (nodes_.empty())
    return;
  Position inv_u {1.0 / u.x, 1.0 / u.y, 1.0 / u.z};
  array<int32_t, MAX_DEPTH> stack;
  array<double, MAX_DEPTH> stack_t;
  int n = 0;
  int32_t i = 0;
  double t = enter(nodes_[0].box, r, inv_u, t_max);
  while (true) {
    if (t <= t_max) {
      const Node& node = nodes_[i];
      if (node.count == 0) {
        // Visit the nearer child first
        int32_t a = i + 1;
        int32_t b = node.first;
        double ta = enter(nodes_[a].box, r, inv_u, t_max);
        double tb = enter(nodes_[b].box, r, inv_u, t_max);
        if (tb < ta) {
          std::swap(a, b);
          std::swap(ta, tb);
        }
        if (tb <= t_max) {
          stack[n] = b;
          stack_t[n++] = tb;
        }
        i = a;
        t = ta;
        continue;
      }
      for (int32_t k = node.first; k < node.first + node.count; ++k) {
        double tk = enter(item_boxes_[k], r, inv_u, t_max);
        if (tk <= t_max)
          f(items_[k], tk);
      }
    }
    if (n == 0)
      return;
    --n;
    i = stack[n];
    t = stack_t[n];
  }
}

} // namespace openmc

#endif // OPENMC_BOX_TREE_H
