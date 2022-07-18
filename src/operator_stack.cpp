#include "openmc/operator_stack.h"

namespace openmc {

/*!
 * Whether the stack has any pushed values
 */
bool OperatorStack::empty() const
{
  return size_ == 0;
}

/*!
 * Push an operator onto the stack
 */
void OperatorStack::push_back(value_type op)
{
  OperatorStack::size_type bit_op = convert_to_bit_notation(op);
  data_ = shl(data_) | bit_op;
  ++size_;
}

/*!
 * Pop a value off the stack
 */
OperatorStack::value_type OperatorStack::pop()
{
  OperatorStack::size_type bit_op = top(data_);
  data_ = shr(data_);
  --size_;
  return convert_to_std_notation(bit_op);
}

/*!
 * Convert operator to bit notation
 * OP_UNION: 1 (0001)
 * OP_INTERSECTION: 2 (0010)
 * OP_COMPLEMENT: 3 (0011)
 */
OperatorStack::size_type OperatorStack::convert_to_bit_notation(value_type op)
{
  if (op == OP_UNION) {
    return BIT_OP_UNION;
  } else if (op == OP_INTERSECTION) {
    return BIT_OP_INTERSECTION;
  } else {
    return BIT_OP_COMPLEMENT;
  }
}

/*!
 * Convert operator from bit notation to int32_t
 */
OperatorStack::value_type OperatorStack::convert_to_std_notation(
  size_type bit_op)
{
  if (bit_op == BIT_OP_UNION) {
    return OP_UNION;
  } else if (bit_op == BIT_OP_INTERSECTION) {
    return OP_INTERSECTION;
  } else {
    return OP_COMPLEMENT;
  }
}

} // namespace openmc
