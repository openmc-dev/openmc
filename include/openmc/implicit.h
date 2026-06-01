#ifndef OPENMC_IMPLICIT_H
#define OPENMC_IMPLICIT_H

#include <cmath>
#include <cstdint>
#include <limits> // For numeric_limits
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

#include "pugixml.hpp"

#include "openmc/position.h"

namespace openmc {

//==============================================================================
// Type aliases
//==============================================================================

// Gradient shares the memory layout of Position (x, y, z doubles).
// Using a named alias makes intent explicit at every call site.
using Gradient = Position;

//==============================================================================
// StepCache
//
// Thread-local epoch counter.  Increment step_cache.epoch once at the
// start of each geometry step.  All Cached nodes across the DAG are
// lazily invalidated at zero cost — no traversal needed.
//==============================================================================

struct CacheEntry {
  uint64_t epoch {std::numeric_limits<uint64_t>::max()};
  Position pos {};
  double val {0.0};
  Gradient grad {0.0, 0.0, 0.0};
};

struct StepCache {
  uint64_t epoch {0};
  // Per-thread cache keyed by node address.
  // Grows to at most one entry per Cached node in the model — bounded and tiny.
  std::unordered_map<const void*, CacheEntry> node_cache;
};

extern thread_local StepCache step_cache;

//==============================================================================
// Implicit
//
// Abstract base class for expression DAG nodes representing a smooth
// 3D function f(x, y, z).  Concrete subclasses are defined in the
// openmc::implicit namespace below.
//
// The DAG is constructed in Python (implicit_function.py) and serialised
// to XML.  On the C++ side it is reconstructed via from_xml_element() and
// evaluated during implicit surface geometry queries.
//==============================================================================

class Implicit {
public:
  virtual ~Implicit() = default;

  virtual std::string expression() const = 0;
  //! Evaluate f(r).
  virtual double evaluate(Position r) const = 0;

  //! Gradient df/dx, df/dy, df/dz.
  //! Computed analytically — no finite differences.
  virtual Gradient gradient(Position r) const = 0;

  //! Evaluate f along the ray: f(r + t*u).
  double along_ray(double t, Position r, Direction u) const;

  //! Lipschitz constant of f on the ray segment [t0, t1]:
  //!   |f(r + a*u) - f(r + b*u)| <= L * |a - b|  for all a,b in [t0,t1].
  //! Computed analytically via interval arithmetic propagation.
  virtual double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const = 0;

  //! Tight lower and upper bounds of f on the ray segment [t0, t1].
  //! Returned as {f_min, f_max}.
  //! Computed jointly because many nodes (Mul, Div, Pow) need both bounds
  //! together to propagate the interval correctly.
  virtual std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const = 0;

  // create a xml element out of the node
  virtual pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const = 0;

  // Convenience wrapper — serialises the whole DAG to an XML string
  std::string to_xml_string() const;

  //! Reconstruct an Implicit DAG from an XML element produced by
  //! Python's to_xml_element().
  //!
  //! \param node       The XML element to parse.
  //! \param cache_map  Maps numeric ids to already-constructed shared nodes.
  //!                   Populated when a <to_cache> element is encountered;
  //!                   looked up when a <from_cache> element is encountered.
  //!                   Passed by reference through the recursion.
  static std::shared_ptr<Implicit> from_xml_element(pugi::xml_node node,
    std::unordered_map<int, std::shared_ptr<Implicit>>& cache_map);

  //! Convenience overload for the root call — creates an empty cache_map.
  static std::shared_ptr<Implicit> from_xml_element(pugi::xml_node node);
};

//==============================================================================
// Concrete node types  (openmc::implicit namespace)
//
// All implementations are in implicit.cpp.
// Nodes follow the expression tree structure of the Python counterparts in
// implicit_function.py.  The openmc::implicit sub-namespace avoids name
// collisions with common identifiers (X, Y, Z, Log, Exp, etc.).
//==============================================================================

namespace implicit {

// ============================================================================
// Terminal nodes
// ============================================================================

class X final : public Implicit {
public:
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;
};

class Y final : public Implicit {
public:
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;
};

class Z final : public Implicit {
public:
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;
};

class Constant final : public Implicit {
public:
  explicit Constant(double value) : value_(value) {}
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;

private:
  double value_;
};

// ============================================================================
// Arithmetic nodes
// ============================================================================

class Add final : public Implicit {
public:
  explicit Add(std::shared_ptr<Implicit> f, std::shared_ptr<Implicit> g)
    : f_(std::move(f)), g_(std::move(g))
  {}
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;

private:
  std::shared_ptr<Implicit> f_, g_;
};

class Sub final : public Implicit {
public:
  explicit Sub(std::shared_ptr<Implicit> f, std::shared_ptr<Implicit> g)
    : f_(std::move(f)), g_(std::move(g))
  {}
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;

private:
  std::shared_ptr<Implicit> f_, g_;
};

class Mul final : public Implicit {
public:
  explicit Mul(std::shared_ptr<Implicit> f, std::shared_ptr<Implicit> g)
    : f_(std::move(f)), g_(std::move(g))
  {}
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;

private:
  std::shared_ptr<Implicit> f_, g_;
};

class Div final : public Implicit {
public:
  explicit Div(std::shared_ptr<Implicit> f, std::shared_ptr<Implicit> g)
    : f_(std::move(f)), g_(std::move(g))
  {}
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;

private:
  std::shared_ptr<Implicit> f_, g_;
};

class Scale final : public Implicit {
public:
  explicit Scale(std::shared_ptr<Implicit> f, double k)
    : f_(std::move(f)), k_(k)
  {}
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;

private:
  std::shared_ptr<Implicit> f_;
  double k_;
};

class Neg final : public Implicit {
public:
  explicit Neg(std::shared_ptr<Implicit> f) : f_(std::move(f)) {}
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;

private:
  std::shared_ptr<Implicit> f_;
};

class Pow final : public Implicit {
public:
  explicit Pow(std::shared_ptr<Implicit> f, int exp)
    : f_(std::move(f)), exp_(exp)
  {
    if (exp_ <= 0)
      throw std::domain_error(
        "Pow: exponent must be a strictly positive integer, got " +
        std::to_string(exp_));
  }
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;

private:
  std::shared_ptr<Implicit> f_;
  int exp_;
};

// ============================================================================
// Transcendental nodes
// ============================================================================

class Sin final : public Implicit {
public:
  explicit Sin(std::shared_ptr<Implicit> arg) : arg_(std::move(arg)) {}
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;

private:
  std::shared_ptr<Implicit> arg_;
};

class Cos final : public Implicit {
public:
  explicit Cos(std::shared_ptr<Implicit> arg) : arg_(std::move(arg)) {}
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;

private:
  std::shared_ptr<Implicit> arg_;
};

class Sqrt final : public Implicit {
public:
  explicit Sqrt(std::shared_ptr<Implicit> arg) : arg_(std::move(arg)) {}
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;

private:
  std::shared_ptr<Implicit> arg_;
};

class Exp final : public Implicit {
public:
  explicit Exp(std::shared_ptr<Implicit> arg) : arg_(std::move(arg)) {}
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;

private:
  std::shared_ptr<Implicit> arg_;
};

class Log final : public Implicit {
public:
  explicit Log(std::shared_ptr<Implicit> arg) : arg_(std::move(arg)) {}
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;

private:
  std::shared_ptr<Implicit> arg_;
};

class Abs final : public Implicit {
public:
  explicit Abs(std::shared_ptr<Implicit> arg) : arg_(std::move(arg)) {}
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;

private:
  std::shared_ptr<Implicit> arg_;
};

// ============================================================================
// Cached node
//
// Memoises evaluate() and gradient() for the duration of one geometry step
// AND one position.  A cache hit requires BOTH a matching epoch AND a
// matching position.
//
// Epoch-only invalidation is INCORRECT: the ray-tracing solver evaluates f
// at many different positions within a single step, so the position check
// is essential for correctness.
//
// compute_lipschitz and compute_f_min_max are NOT cached — they operate on
// intervals, not single points, so the per-point cache does not apply.
// Both methods delegate directly to the child.
//
// Thread safety: Cached has no mutable members. All cache data lives in
// thread_local step_cache.node_cache, keyed by node address. Each thread
// maintains a fully independent cache — no shared mutable state, no data race.
//
// The map grows to at most one entry per Cached node in the model and is
// never cleared — entries are simply overwritten when the epoch or position
// changes.
//
// WARNING: the SAME shared_ptr instance must be used wherever this
// subexpression appears in the DAG.  Constructing two separate Cached nodes
// wrapping the same expression gives two independent caches and no XML
// sharing.  See the Python-side docstring in implicit_function.py for the
// full pitfall and correct usage pattern.
// ============================================================================

class Cached final : public Implicit {
public:
  explicit Cached(std::shared_ptr<Implicit> child) : child_(std::move(child)) {}
  std::string expression() const override;
  double evaluate(Position r) const override;
  Gradient gradient(Position r) const override;
  double compute_lipschitz(
    Position r, Direction u, double t0, double t1) const override;
  std::pair<double, double> compute_f_min_max(
    Position r, Direction u, double t0, double t1) const override;
  pugi::xml_node to_xml_element(pugi::xml_node parent,
    std::unordered_map<const Implicit*, int>& cache_map) const override;

private:
  //! Recompute value and gradient if epoch or position changed.
  void refresh(Position r) const;

  std::shared_ptr<Implicit> child_;
};

} // namespace implicit
} // namespace openmc
#endif // OPENMC_IMPLICIT_H
