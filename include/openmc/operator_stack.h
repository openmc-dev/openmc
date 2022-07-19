#ifndef OPENMC_OPERATOR_STACK
#define OPENMC_OPERATOR_STACK

#include <cstdint>

#include "cell.h"

namespace openmc {

//==============================================================================

class OperatorStack {
public:
  //----------------------------------------------------------------------------
  // Typedefs
  using value_type = int32_t;
  using size_type = uint32_t;

  //----------------------------------------------------------------------------
  // Constructors, destructors

  //! Default constructor
  OperatorStack() {}

  //----------------------------------------------------------------------------
  // Accessors

  //! Greatest number of operators allowed on the stack, each operator takes two
  //! bits
  static constexpr size_type max_stack_depth() { return sizeof(size_type) * 4; }

  // Whether any elements exist
  bool empty() const { return data_ == size_type(0); }

  //----------------------------------------------------------------------------
  // Mutators

  // Push an operator to the stack
  void push_back(value_type op) { data_ = shl(data_) | twobit_operator(op); }

  // Pop a value off the stack
  value_type pop()
  {
    auto bit_op = top(data_);
    data_ = shr(data_);
    return signed_operator(bit_op);
  }

private:
  //----------------------------------------------------------------------------
  // Private methods
  size_type twobit_operator(value_type op) const { return op - OP_UNION + 1; }

  value_type signed_operator(size_type bit_op) const
  {
    return bit_op + OP_UNION - 1;
  }

  //! Get the top two bits
  static constexpr size_type top(size_type val) { return val & 0b11u; }

  //! Shift right by two
  static constexpr size_type shr(size_type val) { return val >> 2; }

  //! Shift left by two
  static constexpr size_type shl(size_type val) { return val << 2; }

  //----------------------------------------------------------------------------
  // Data
  size_type data_ {0};
};

} // namespace openmc

#endif // OPENMC_OPERATOR_STACK
