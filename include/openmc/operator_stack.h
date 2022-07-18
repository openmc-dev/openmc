#ifndef OPENMC_OPERATOR_STACK
#define OPENMC_OPERATOR_STACK

#include <cstdint>

#include "cell.h"

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

enum BitOperators : uint16_t {
  BIT_OP_UNION = 1,
  BIT_OP_INTERSECTION = 2,
  BIT_OP_COMPLEMENT = 3
};

//==============================================================================
//==============================================================================

class OperatorStack {
public:
  //@{
  //! Typedefs
  using value_type = int32_t;
  using size_type = uint16_t;
  //}

  //! Default constructor
  OperatorStack() {}

  //! Greatest number of operators allowed on the stack, each operator takes two
  //! bits
  static constexpr size_type max_stack_depth() { return sizeof(size_type) * 4; }

  //// ACCESSORS ////

  //! Number of elements on the stack
  size_type size() const { return size_; }

  // Whether any elements exist
  bool empty() const;

  //// MUTATORS ////

  // Push an operator to the stack
  void push_back(value_type op);

  // Pop a value off the stack
  value_type pop();

  // Convert a operator to bit operator
  size_type convert_to_bit_notation(value_type op);

  // Convert a bit operator to operator
  value_type convert_to_std_notation(size_type bit_op);

private:
  //// DATA ////

  size_type data_ {0};
  size_type size_ {0};

  //! Get the top two bits
  static constexpr size_type top(size_type val) { return val & 3; }

  //! Shift right by two
  static constexpr size_type shr(size_type val) { return val >> 2; }

  //! Shift left by two
  static constexpr size_type shl(size_type val) { return val << 2; }
};

} // namespace openmc

#endif // OPENMC_OPERATOR_STACK
