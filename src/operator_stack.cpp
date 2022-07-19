#include "openmc/operator_stack.h"

namespace openmc {

/*!
 * Whether the stack has any pushed values
 */
bool OperatorStack::empty() const
{
  return data_ == 0;
}

/*!
 * Push an operator onto the stack
 */
void OperatorStack::push_back(value_type op)
{
  data_ = shl(data_) | convert_to_bit_notation(op);
}

/*!
 * Pop a value off the stack
 */
OperatorStack::value_type OperatorStack::pop()
{
  auto bit_op = top(data_);
  data_ = shr(data_);
  return convert_to_std_notation(bit_op);
}

/*!
 * Convert operator to bit notation
 * OP_UNION: 1 (0001)
 * OP_INTERSECTION: 2 (0010)
 * OP_COMPLEMENT: 3 (0011)
 */
OperatorStack::size_type OperatorStack::convert_to_bit_notation(
  value_type op) const
{
  switch (op) {
  case OP_UNION:
    return BIT_OP_UNION;
  case OP_INTERSECTION:
    return BIT_OP_INTERSECTION;
  case OP_COMPLEMENT:
    return BIT_OP_COMPLEMENT;
  default:
    return 0;
  }
}

/*!
 * Convert operator from bit notation to int32_t
 */
OperatorStack::value_type OperatorStack::convert_to_std_notation(
  size_type bit_op) const
{
  switch (bit_op) {
  case BIT_OP_UNION:
    return OP_UNION;
  case BIT_OP_INTERSECTION:
    return OP_INTERSECTION;
  case BIT_OP_COMPLEMENT:
    return OP_COMPLEMENT;
  default:
    return 0;
  }
}

} // namespace openmc
