#pragma once

#ifdef __CUDACC__
#include "memory.h"

template<class Key, class Value, class Allocator = UnifiedAllocator>
class unordered_map {};

#else

#include <unordered_map>
template<class Key, class Value>
using unordered_map

#endif
