#include "openmc/implicit.h"

#include <algorithm> // std::min, std::max
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace openmc {

//==============================================================================
// Cache
//==============================================================================

thread_local StepCache step_cache;

//==============================================================================
// Implicit
//==============================================================================

double Implicit::along_ray(double t, Position r, Direction u) const
{
  return evaluate({r.x + t * u.x, r.y + t * u.y, r.z + t * u.z});
}

std::shared_ptr<Implicit> Implicit::from_xml_element(pugi::xml_node node,
  std::unordered_map<int, std::shared_ptr<Implicit>>& cache_map)
{
  std::string tag = node.name();

  // from_cache: look up a previously registered shared node.
  // Handled before child parsing - there are no children to recurse into.
  if (tag == "from_cache") {
    int id = node.attribute("id").as_int();
    auto it = cache_map.find(id);
    if (it == cache_map.end())
      throw std::runtime_error(
        "Implicit::from_xml_element: from_cache id=" + std::to_string(id) +
        " has no matching to_cache.");
    return it->second;
  }

  // Recursively parse all element children.
  std::vector<std::shared_ptr<Implicit>> children;
  for (auto child : node.children()) {
    if (child.type() == pugi::node_element)
      children.push_back(from_xml_element(child, cache_map));
  }

  // to_cache: register the single child node by id, return a Cached wrapper.
  if (tag == "to_cache") {
    int id = node.attribute("id").as_int();
    auto cached = std::make_shared<implicit::Cached>(children[0]);
    cache_map[id] = cached;
    return cached;
  }

  // Dispatch table - mirrors Python's from_xml_element.
  if (tag == "x")
    return std::make_shared<implicit::X>();
  if (tag == "y")
    return std::make_shared<implicit::Y>();
  if (tag == "z")
    return std::make_shared<implicit::Z>();
  if (tag == "constant") {
    double v = node.attribute("value").as_double();
    return std::make_shared<implicit::Constant>(v);
  }
  if (tag == "add")
    return std::make_shared<implicit::Add>(children[0], children[1]);
  if (tag == "sub")
    return std::make_shared<implicit::Sub>(children[0], children[1]);
  if (tag == "mul")
    return std::make_shared<implicit::Mul>(children[0], children[1]);
  if (tag == "div")
    return std::make_shared<implicit::Div>(children[0], children[1]);
  if (tag == "scale") {
    double k = node.attribute("value").as_double();
    return std::make_shared<implicit::Scale>(children[0], k);
  }
  if (tag == "neg")
    return std::make_shared<implicit::Neg>(children[0]);
  if (tag == "pow") {
    int exp = node.attribute("value").as_int();
    return std::make_shared<implicit::Pow>(children[0], exp);
  }
  if (tag == "sin")
    return std::make_shared<implicit::Sin>(children[0]);
  if (tag == "cos")
    return std::make_shared<implicit::Cos>(children[0]);
  if (tag == "sqrt")
    return std::make_shared<implicit::Sqrt>(children[0]);
  if (tag == "exp")
    return std::make_shared<implicit::Exp>(children[0]);
  if (tag == "log")
    return std::make_shared<implicit::Log>(children[0]);
  if (tag == "abs")
    return std::make_shared<implicit::Abs>(children[0]);

  throw std::runtime_error(
    "Implicit::from_xml_element: unknown tag '" + tag + "'.");
}

std::shared_ptr<Implicit> Implicit::from_xml_element(pugi::xml_node node)
{
  std::unordered_map<int, std::shared_ptr<Implicit>> cache_map;
  return from_xml_element(node, cache_map);
}

std::string Implicit::to_xml_string() const
{
  pugi::xml_document doc;
  auto func_node = doc.append_child("function");
  std::unordered_map<const Implicit*, int> cache_map;
  to_xml_element(func_node, cache_map);
  std::ostringstream oss;
  doc.save(oss);
  return oss.str();
}

//==============================================================================
// Concrete node types  (openmc::implicit namespace)
//==============================================================================

namespace implicit {

//==============================================================================
// Helper functions for the product of two intervals [a,b] * [c,d]
//==============================================================================

static std::pair<double, double> interval_mul(
  double a, double b, double c, double d)
{
  double p1 = a * c, p2 = a * d, p3 = b * c, p4 = b * d;
  return {std::min({p1, p2, p3, p4}), std::max({p1, p2, p3, p4})};
}

//==============================================================================
// X
//==============================================================================

std::string X::expression() const
{
  return "X";
}
double X::evaluate(Position r) const
{
  return r.x;
}
Gradient X::gradient(Position r) const
{
  return {1.0, 0.0, 0.0};
}
double X::compute_lipschitz(Position r, Direction u, double t0, double t1) const
{
  return std::abs(u.x);
}
std::pair<double, double> X::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  double f0 = r.x + t0 * u.x;
  double f1 = r.x + t1 * u.x;
  return {std::min(f0, f1), std::max(f0, f1)};
}
pugi::xml_node X::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("x");
  return node;
}

//==============================================================================
// Y
//==============================================================================

std::string Y::expression() const
{
  return "Y";
}
double Y::evaluate(Position r) const
{
  return r.y;
}
Gradient Y::gradient(Position r) const
{
  return {0.0, 1.0, 0.0};
}
double Y::compute_lipschitz(Position r, Direction u, double t0, double t1) const
{
  return std::abs(u.y);
}
std::pair<double, double> Y::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  double f0 = r.y + t0 * u.y;
  double f1 = r.y + t1 * u.y;
  return {std::min(f0, f1), std::max(f0, f1)};
}
pugi::xml_node Y::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("y");
  return node;
}

//==============================================================================
// Z
//==============================================================================

std::string Z::expression() const
{
  return "Z";
}
double Z::evaluate(Position r) const
{
  return r.z;
}
Gradient Z::gradient(Position r) const
{
  return {0.0, 0.0, 1.0};
}
double Z::compute_lipschitz(Position r, Direction u, double t0, double t1) const
{
  return std::abs(u.z);
}
std::pair<double, double> Z::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  double f0 = r.z + t0 * u.z;
  double f1 = r.z + t1 * u.z;
  return {std::min(f0, f1), std::max(f0, f1)};
}
pugi::xml_node Z::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("z");
  return node;
}

//==============================================================================
// Constant
//==============================================================================

std::string Constant::expression() const
{
  return std::to_string(value_);
}
double Constant::evaluate(Position r) const
{
  return value_;
}
Gradient Constant::gradient(Position r) const
{
  return {0.0, 0.0, 0.0};
}
double Constant::compute_lipschitz(
  Position r, Direction u, double t0, double t1) const
{
  return 0.;
}
std::pair<double, double> Constant::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  return {value_, value_};
}
pugi::xml_node Constant::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("constant");
  node.append_attribute("value") = value_;
  return node;
}

//==============================================================================
// Add
//==============================================================================

std::string Add::expression() const
{
  return f_->expression() + " + " + g_->expression();
}
double Add::evaluate(Position r) const
{
  return f_->evaluate(r) + g_->evaluate(r);
}
Gradient Add::gradient(Position r) const
{
  return f_->gradient(r) + g_->gradient(r);
}
double Add::compute_lipschitz(
  Position r, Direction u, double t0, double t1) const
{
  return f_->compute_lipschitz(r, u, t0, t1) +
         g_->compute_lipschitz(r, u, t0, t1);
}
std::pair<double, double> Add::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  auto [f_min, f_max] = f_->compute_f_min_max(r, u, t0, t1);
  auto [g_min, g_max] = g_->compute_f_min_max(r, u, t0, t1);
  return {f_min + g_min, f_max + g_max};
}
pugi::xml_node Add::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("add");
  f_->to_xml_element(node, cache_map);
  g_->to_xml_element(node, cache_map);
  return node;
}

//==============================================================================
// Sub
//==============================================================================

std::string Sub::expression() const
{
  return f_->expression() + " - " + g_->expression();
}
double Sub::evaluate(Position r) const
{
  return f_->evaluate(r) - g_->evaluate(r);
}
Gradient Sub::gradient(Position r) const
{
  return f_->gradient(r) - g_->gradient(r);
}
double Sub::compute_lipschitz(
  Position r, Direction u, double t0, double t1) const
{
  return f_->compute_lipschitz(r, u, t0, t1) +
         g_->compute_lipschitz(r, u, t0, t1);
}
std::pair<double, double> Sub::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  auto [f_min, f_max] = f_->compute_f_min_max(r, u, t0, t1);
  auto [g_min, g_max] = g_->compute_f_min_max(r, u, t0, t1);
  return {f_min - g_max, f_max - g_min};
}
pugi::xml_node Sub::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("sub");
  f_->to_xml_element(node, cache_map);
  g_->to_xml_element(node, cache_map);
  return node;
}

//==============================================================================
// Mul
//==============================================================================

std::string Mul::expression() const
{
  return f_->expression() + " * " + g_->expression();
}
double Mul::evaluate(Position r) const
{
  return f_->evaluate(r) * g_->evaluate(r);
}
Gradient Mul::gradient(Position r) const
{
  double fv = f_->evaluate(r);
  double gv = g_->evaluate(r);
  return gv * f_->gradient(r) + fv * g_->gradient(r);
}
double Mul::compute_lipschitz(
  Position r, Direction u, double t0, double t1) const
{
  double lipschitz_f = f_->compute_lipschitz(r, u, t0, t1);
  double lipschitz_g = g_->compute_lipschitz(r, u, t0, t1);
  auto [f_min, f_max] = f_->compute_f_min_max(r, u, t0, t1);
  auto [g_min, g_max] = g_->compute_f_min_max(r, u, t0, t1);
  double max_abs_f = std::max(std::abs(f_min), std::abs(f_max));
  double max_abs_g = std::max(std::abs(g_min), std::abs(g_max));
  return max_abs_f * lipschitz_g + max_abs_g * lipschitz_f;
}
std::pair<double, double> Mul::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  auto [f_min, f_max] = f_->compute_f_min_max(r, u, t0, t1);
  auto [g_min, g_max] = g_->compute_f_min_max(r, u, t0, t1);
  auto [p_min, p_max] = interval_mul(f_min, f_max, g_min, g_max);
  return {p_min, p_max};
}
pugi::xml_node Mul::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("mul");
  f_->to_xml_element(node, cache_map);
  g_->to_xml_element(node, cache_map);
  return node;
}

//==============================================================================
// Div
//==============================================================================

std::string Div::expression() const
{
  return f_->expression() + " / " + g_->expression();
}
double Div::evaluate(Position r) const
{
  double denom = g_->evaluate(r);
  if (denom == 0.0)
    throw std::domain_error("Implicit Div: zero denominator at r=(" +
                            std::to_string(r.x) + ", " + std::to_string(r.y) +
                            ", " + std::to_string(r.z) + ") in expression " +
                            g_->expression());
  return f_->evaluate(r) / denom;
}
Gradient Div::gradient(Position r) const
{
  double fv = f_->evaluate(r);
  double gv = g_->evaluate(r);
  Gradient gf = f_->gradient(r);
  Gradient gg = g_->gradient(r);
  if (gv == 0.0) {
    throw std::domain_error("Div::gradient: argument equals zero at r=(" +
                            std::to_string(r.x) + ", " + std::to_string(r.y) +
                            ", " + std::to_string(r.z) + ") in expression " +
                            g_->expression());
  }
  double inv_g2 = 1.0 / (gv * gv);
  return (gf * gv - fv * gg) * inv_g2;
}
double Div::compute_lipschitz(
  Position r, Direction u, double t0, double t1) const
{
  double lipschitz_f = f_->compute_lipschitz(r, u, t0, t1);
  double lipschitz_g = g_->compute_lipschitz(r, u, t0, t1);
  auto [f_min, f_max] = f_->compute_f_min_max(r, u, t0, t1);
  auto [g_min, g_max] = g_->compute_f_min_max(r, u, t0, t1);
  double min_abs_g = std::min(std::abs(g_min), std::abs(g_max));
  double max_abs_f = std::max(std::abs(f_min), std::abs(f_max));
  double max_abs_g = std::max(std::abs(g_min), std::abs(g_max));
  if (min_abs_g > 0.0) {
    return (lipschitz_f * max_abs_g + max_abs_f * lipschitz_g) /
           (min_abs_g * min_abs_g);
  } else {
    throw std::domain_error(
      "Implicit Div: zero denominator on ray between r=(" +
      std::to_string(r.x + t0 * u.x) + ", " + std::to_string(r.y + t0 * u.y) +
      ", " + std::to_string(r.z + t0 * u.z) + ") and r=(" +
      std::to_string(r.x + t1 * u.x) + ", " + std::to_string(r.y + t1 * u.y) +
      ", " + std::to_string(r.z + t1 * u.z) + ") in expression " +
      g_->expression());
  }
}
std::pair<double, double> Div::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  auto [f_min, f_max] = f_->compute_f_min_max(r, u, t0, t1);
  auto [g_min, g_max] = g_->compute_f_min_max(r, u, t0, t1);
  double inv_min, inv_max;
  if (g_min > 0.0) {
    // g does not contain zero - compute 1/[g_min, g_max]
    inv_min = 1.0 / g_max;
    inv_max = 1.0 / g_min;
  } else if (g_max < 0.0) {
    // g does not contain zero - compute 1/[g_min, g_max]
    inv_min = 1.0 / g_min;
    inv_max = 1.0 / g_max;
  } else {
    // g straddles zero - not a lipschitz function.
    throw std::domain_error(
      "Implicit Div: zero denominator on ray between r=(" +
      std::to_string(r.x + t0 * u.x) + ", " + std::to_string(r.y + t0 * u.y) +
      ", " + std::to_string(r.z + t0 * u.z) + ") and r=(" +
      std::to_string(r.x + t1 * u.x) + ", " + std::to_string(r.y + t1 * u.y) +
      ", " + std::to_string(r.z + t1 * u.z) + ") in expression " +
      g_->expression());
  }
  auto [p_min, p_max] = interval_mul(f_min, f_max, inv_min, inv_max);
  return {p_min, p_max};
}
pugi::xml_node Div::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("div");
  f_->to_xml_element(node, cache_map);
  g_->to_xml_element(node, cache_map);
  return node;
}

//==============================================================================
// Scale
//==============================================================================

std::string Scale::expression() const
{
  return std::to_string(k_) + " * " + f_->expression();
}
double Scale::evaluate(Position r) const
{
  return k_ * f_->evaluate(r);
}
Gradient Scale::gradient(Position r) const
{
  return k_ * f_->gradient(r);
}
double Scale::compute_lipschitz(
  Position r, Direction u, double t0, double t1) const
{
  return std::abs(k_) * f_->compute_lipschitz(r, u, t0, t1);
}
std::pair<double, double> Scale::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  auto [f_min, f_max] = f_->compute_f_min_max(r, u, t0, t1);
  if (k_ >= 0.0) {
    return {k_ * f_min, k_ * f_max};
  } else {
    return {k_ * f_max, k_ * f_min};
  }
}
pugi::xml_node Scale::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("scale");
  node.append_attribute("value") = k_;
  f_->to_xml_element(node, cache_map);
  return node;
}

//==============================================================================
// Neg
//==============================================================================

std::string Neg::expression() const
{
  return "-" + f_->expression();
}
double Neg::evaluate(Position r) const
{
  return -f_->evaluate(r);
}
Gradient Neg::gradient(Position r) const
{
  return -f_->gradient(r);
}
double Neg::compute_lipschitz(
  Position r, Direction u, double t0, double t1) const
{
  return f_->compute_lipschitz(r, u, t0, t1);
}
std::pair<double, double> Neg::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  auto [f_min, f_max] = f_->compute_f_min_max(r, u, t0, t1);
  return {-f_max, -f_min};
}
pugi::xml_node Neg::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("neg");
  f_->to_xml_element(node, cache_map);
  return node;
}

//==============================================================================
// Pow
//==============================================================================

std::string Pow::expression() const
{
  return f_->expression() + " ** " + std::to_string(exp_);
}
double Pow::evaluate(Position r) const
{
  return std::pow(f_->evaluate(r), exp_);
}
Gradient Pow::gradient(Position r) const
{
  // Chain rule: (f^n)' = n * f^(n-1) * f'
  double fv = f_->evaluate(r);
  Gradient gf = f_->gradient(r);
  double coeff = exp_ * std::pow(fv, exp_ - 1.0);
  return coeff * gf;
}
double Pow::compute_lipschitz(
  Position r, Direction u, double t0, double t1) const
{
  double lipschitz_f = f_->compute_lipschitz(r, u, t0, t1);
  auto [f_min, f_max] = f_->compute_f_min_max(r, u, t0, t1);
  double max_abs_f = std::max(std::abs(f_min), std::abs(f_max));
  return std::abs(exp_) * std::pow(max_abs_f, exp_ - 1.0) * lipschitz_f;
}
std::pair<double, double> Pow::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  auto [f_min, f_max] = f_->compute_f_min_max(r, u, t0, t1);
  if (exp_ % 2 != 0) {
    // Odd exponent: monotone increasing everywhere
    return {std::pow(f_min, exp_), std::pow(f_max, exp_)};
  } else {
    // Even exponent: minimum at 0 if interval straddles 0
    double p_min = std::pow(f_min, exp_);
    double p_max = std::pow(f_max, exp_);
    if (f_min <= 0.0 && f_max >= 0.0)
      return {0.0, std::max(p_min, p_max)};
    else
      return {std::min(p_min, p_max), std::max(p_min, p_max)};
  }
}
pugi::xml_node Pow::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("pow");
  node.append_attribute("value") = exp_;
  f_->to_xml_element(node, cache_map);
  return node;
}

//==============================================================================
// Sin
//==============================================================================

std::string Sin::expression() const
{
  return "Sin(" + arg_->expression() + ")";
}
double Sin::evaluate(Position r) const
{
  return std::sin(arg_->evaluate(r));
}
Gradient Sin::gradient(Position r) const
{
  double c = std::cos(arg_->evaluate(r));
  Gradient g = arg_->gradient(r);
  return c * g;
}
double Sin::compute_lipschitz(
  Position r, Direction u, double t0, double t1) const
{
  return arg_->compute_lipschitz(r, u, t0, t1);
}
std::pair<double, double> Sin::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  auto [arg_min, arg_max] = arg_->compute_f_min_max(r, u, t0, t1);
  double sin_arg_min = std::sin(arg_min);
  double sin_arg_max = std::sin(arg_max);
  double f_min = std::min(sin_arg_min, sin_arg_max);
  double f_max = std::max(sin_arg_min, sin_arg_max);

  // sin reaches +1 at π/2 + 2kπ
  double period_arg_min = std::ceil((arg_min - M_PI / 2.0) / (2.0 * M_PI));
  double period_arg_max = std::floor((arg_max - M_PI / 2.0) / (2.0 * M_PI));
  if (period_arg_min <= period_arg_max)
    f_max = 1.0;

  // sin reaches -1 at -π/2 + 2kπ
  period_arg_min = std::ceil((arg_min + M_PI / 2.0) / (2.0 * M_PI));
  period_arg_max = std::floor((arg_max + M_PI / 2.0) / (2.0 * M_PI));
  if (period_arg_min <= period_arg_max)
    f_min = -1.0;

  return {f_min, f_max};
}
pugi::xml_node Sin::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("sin");
  arg_->to_xml_element(node, cache_map);
  return node;
}

//==============================================================================
// Cos
//==============================================================================

std::string Cos::expression() const
{
  return "Cos(" + arg_->expression() + ")";
}
double Cos::evaluate(Position r) const
{
  return std::cos(arg_->evaluate(r));
}
Gradient Cos::gradient(Position r) const
{
  double c = -std::sin(arg_->evaluate(r));
  Gradient g = arg_->gradient(r);
  return c * g;
}
double Cos::compute_lipschitz(
  Position r, Direction u, double t0, double t1) const
{
  return arg_->compute_lipschitz(r, u, t0, t1);
}
std::pair<double, double> Cos::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  auto [arg_min, arg_max] = arg_->compute_f_min_max(r, u, t0, t1);
  double cos_arg_min = std::cos(arg_min);
  double cos_arg_max = std::cos(arg_max);
  double f_min = std::min(cos_arg_min, cos_arg_max);
  double f_max = std::max(cos_arg_min, cos_arg_max);

  // cos reaches +1 at 2kπ
  double period_arg_min = std::ceil(arg_min / (2.0 * M_PI));
  double period_arg_max = std::floor(arg_max / (2.0 * M_PI));
  if (period_arg_min <= period_arg_max)
    f_max = 1.0;

  // cos reaches -1 at π + 2kπ
  period_arg_min = std::ceil((arg_min - M_PI) / (2.0 * M_PI));
  period_arg_max = std::floor((arg_max - M_PI) / (2.0 * M_PI));
  if (period_arg_min <= period_arg_max)
    f_min = -1.0;

  return {f_min, f_max};
}
pugi::xml_node Cos::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("cos");
  arg_->to_xml_element(node, cache_map);
  return node;
}

//==============================================================================
// Sqrt
//==============================================================================

std::string Sqrt::expression() const
{
  return "Sqrt(" + arg_->expression() + ")";
}
double Sqrt::evaluate(Position r) const
{
  return std::sqrt(arg_->evaluate(r));
}
Gradient Sqrt::gradient(Position r) const
{
  double argv = arg_->evaluate(r);
  if (argv <= 0.0) {
    throw std::domain_error(
      "Sqrt::gradient: argument equals zero or negative r=(" +
      std::to_string(r.x) + ", " + std::to_string(r.y) + ", " +
      std::to_string(r.z) + ") in expression " + arg_->expression());
  }
  double c = 0.5 / std::sqrt(argv);
  Gradient g = arg_->gradient(r);
  return c * g;
}
double Sqrt::compute_lipschitz(
  Position r, Direction u, double t0, double t1) const
{
  // L(√f) = L(f) / (2 * √f_min) — derivative 1/(2√f) is largest at f_min
  double lipschitz_arg = arg_->compute_lipschitz(r, u, t0, t1);
  auto [arg_min, arg_max] = arg_->compute_f_min_max(r, u, t0, t1);
  if (arg_min <= 0.0)
    throw std::domain_error(
      "Sqrt::compute_lipschitz: argument reaches zero or negative on ray "
      "between r=(" +
      std::to_string(r.x + t0 * u.x) + ", " + std::to_string(r.y + t0 * u.y) +
      ", " + std::to_string(r.z + t0 * u.z) + ") and r=(" +
      std::to_string(r.x + t1 * u.x) + ", " + std::to_string(r.y + t1 * u.y) +
      ", " + std::to_string(r.z + t1 * u.z) + ") in expression " +
      arg_->expression());
  return lipschitz_arg / (2.0 * std::sqrt(arg_min));
}
std::pair<double, double> Sqrt::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  auto [arg_min, arg_max] = arg_->compute_f_min_max(r, u, t0, t1);
  if (arg_min < 0.0)
    throw std::domain_error(
      "Sqrt::compute_f_min_max: argument is negative on ray between r=(" +
      std::to_string(r.x + t0 * u.x) + ", " + std::to_string(r.y + t0 * u.y) +
      ", " + std::to_string(r.z + t0 * u.z) + ") and r=(" +
      std::to_string(r.x + t1 * u.x) + ", " + std::to_string(r.y + t1 * u.y) +
      ", " + std::to_string(r.z + t1 * u.z) + ") in expression " +
      arg_->expression());
  return {std::sqrt(arg_min), std::sqrt(arg_max)};
}
pugi::xml_node Sqrt::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("sqrt");
  arg_->to_xml_element(node, cache_map);
  return node;
}

//==============================================================================
// Exp
//==============================================================================

std::string Exp::expression() const
{
  return "Exp(" + arg_->expression() + ")";
}
double Exp::evaluate(Position r) const
{
  return std::exp(arg_->evaluate(r));
}
Gradient Exp::gradient(Position r) const
{
  double c = std::exp(arg_->evaluate(r));
  Gradient g = arg_->gradient(r);
  return c * g;
}
double Exp::compute_lipschitz(
  Position r, Direction u, double t0, double t1) const
{
  double lipschitz_arg = arg_->compute_lipschitz(r, u, t0, t1);
  auto [arg_min, arg_max] = arg_->compute_f_min_max(r, u, t0, t1);
  return lipschitz_arg * std::exp(arg_max);
}
std::pair<double, double> Exp::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  auto [arg_min, arg_max] = arg_->compute_f_min_max(r, u, t0, t1);
  return {std::exp(arg_min), std::exp(arg_max)};
}
pugi::xml_node Exp::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("exp");
  arg_->to_xml_element(node, cache_map);
  return node;
}

//==============================================================================
// Log
//==============================================================================

std::string Log::expression() const
{
  return "Log(" + arg_->expression() + ")";
}
double Log::evaluate(Position r) const
{
  return std::log(arg_->evaluate(r));
}
Gradient Log::gradient(Position r) const
{
  double argv = arg_->evaluate(r);
  if (argv <= 0.0) {
    throw std::domain_error(
      "Log::gradient: argument equals zero or negative r=(" +
      std::to_string(r.x) + ", " + std::to_string(r.y) + ", " +
      std::to_string(r.z) + ") in expression " + arg_->expression());
  }
  double c = 1 / argv;
  Gradient g = arg_->gradient(r);
  return c * g;
}
double Log::compute_lipschitz(
  Position r, Direction u, double t0, double t1) const
{
  double lipschitz_arg = arg_->compute_lipschitz(r, u, t0, t1);
  auto [arg_min, arg_max] = arg_->compute_f_min_max(r, u, t0, t1);
  if (arg_min <= 0.0)
    throw std::domain_error(
      "Log::compute_lipschitz: argument reaches zero or negative on ray "
      "between r=(" +
      std::to_string(r.x + t0 * u.x) + ", " + std::to_string(r.y + t0 * u.y) +
      ", " + std::to_string(r.z + t0 * u.z) + ") and r=(" +
      std::to_string(r.x + t1 * u.x) + ", " + std::to_string(r.y + t1 * u.y) +
      ", " + std::to_string(r.z + t1 * u.z) + ") in expression " +
      arg_->expression());
  return lipschitz_arg / arg_min;
}
std::pair<double, double> Log::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  auto [arg_min, arg_max] = arg_->compute_f_min_max(r, u, t0, t1);
  if (arg_min <= 0.0)
    throw std::domain_error(
      "Log::compute_f_min_max: argument is zero or negative on ray between "
      "r=(" +
      std::to_string(r.x + t0 * u.x) + ", " + std::to_string(r.y + t0 * u.y) +
      ", " + std::to_string(r.z + t0 * u.z) + ") and r=(" +
      std::to_string(r.x + t1 * u.x) + ", " + std::to_string(r.y + t1 * u.y) +
      ", " + std::to_string(r.z + t1 * u.z) + ") in expression " +
      arg_->expression());
  return {std::log(arg_min), std::log(arg_max)};
}
pugi::xml_node Log::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("log");
  arg_->to_xml_element(node, cache_map);
  return node;
}

//==============================================================================
// Abs
//==============================================================================

std::string Abs::expression() const
{
  return "|" + arg_->expression() + "|";
}
double Abs::evaluate(Position r) const
{
  return std::abs(arg_->evaluate(r));
}
Gradient Abs::gradient(Position r) const
{
  double c = std::copysign(1.0, arg_->evaluate(r));
  Gradient g = arg_->gradient(r);
  return c * g;
}
double Abs::compute_lipschitz(
  Position r, Direction u, double t0, double t1) const
{
  return arg_->compute_lipschitz(r, u, t0, t1);
}
std::pair<double, double> Abs::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  auto [arg_min, arg_max] = arg_->compute_f_min_max(r, u, t0, t1);
  if (arg_min >= 0.0)
    return {arg_min, arg_max}; // entirely non-negative
  else if (arg_max <= 0.0)
    return {-arg_max, -arg_min}; // entirely non-positive
  else
    return {0.0, std::max(-arg_min, arg_max)}; // straddles zero
}
pugi::xml_node Abs::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto node = parent.append_child("abs");
  arg_->to_xml_element(node, cache_map);
  return node;
}

//==============================================================================
// Cached
//==============================================================================

void Cached::refresh(Position r) const
{
  auto& entry = step_cache.node_cache[this]; // creates default entry if absent

  if (entry.epoch != step_cache.epoch || r.x != entry.pos.x ||
      r.y != entry.pos.y || r.z != entry.pos.z) {
    entry.val = child_->evaluate(r);
    entry.grad = child_->gradient(r);
    entry.epoch = step_cache.epoch;
    entry.pos = r;
  }
}
std::string Cached::expression() const
{
  return "@[" + child_->expression() + "]";
}
double Cached::evaluate(Position r) const
{
  refresh(r);
  return step_cache.node_cache.at(this).val;
}
Gradient Cached::gradient(Position r) const
{
  refresh(r);
  return step_cache.node_cache.at(this).grad;
}
double Cached::compute_lipschitz(
  Position r, Direction u, double t0, double t1) const
{
  return child_->compute_lipschitz(r, u, t0, t1);
}
std::pair<double, double> Cached::compute_f_min_max(
  Position r, Direction u, double t0, double t1) const
{
  return child_->compute_f_min_max(r, u, t0, t1);
}
pugi::xml_node Cached::to_xml_element(pugi::xml_node parent,
  std::unordered_map<const Implicit*, int>& cache_map) const
{
  auto it = cache_map.find(this);
  if (it != cache_map.end()) {
    // Already emitted — write a back-reference
    auto node = parent.append_child("from_cache");
    node.append_attribute("id") = it->second;
    return node;
  }
  // First visit — register and emit the full subtree
  int id = static_cast<int>(cache_map.size());
  cache_map[this] = id;
  auto node = parent.append_child("to_cache");
  node.append_attribute("id") = id;
  child_->to_xml_element(node, cache_map);
  return node;
}

} // namespace implicit
} // namespace openmc