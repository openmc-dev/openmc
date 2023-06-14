#include "openmc/aabb.h"
#include <stdint.h>
#include <float.h>
#include <algorithm>

namespace openmc {

inline vec3 vmin(vec3 lhs, vec3 rhs) {
    vec3 res;
    res.x = std::min(lhs.x, rhs.x);
    res.y = std::min(lhs.y, rhs.y);
    res.z = std::min(lhs.z, rhs.z);
    return res;
}

inline vec3 vmax(vec3 lhs, vec3 rhs) {
    vec3 res;
    res.x = std::max(lhs.x, rhs.x);
    res.y = std::max(lhs.y, rhs.y);
    res.z = std::max(lhs.z, rhs.z);
    return res;
}

AABB::AABB() : min(FLT_MAX, FLT_MAX, FLT_MAX), max(-FLT_MAX, -FLT_MAX, -FLT_MAX) {}

AABB::AABB(const vec3& mi, const vec3& ma) : min(mi), max(ma) {}

void AABB::reset() {
	AABB clear;
	*this = clear;
}

void AABB::extend_max(const vec3& val) {
	max = vmax(max, val);
}

void AABB::extend_min(const vec3& val) {
	min = vmin(min, val);
}

void AABB::extend(const vec3& pos) {
	extend_max(pos);
	extend_min(pos);
}

float AABB::surface_area() const {
	vec3 sideLengths = max - min;

	return
		2 * (sideLengths.x * (sideLengths.y + sideLengths.z) +
		sideLengths.y *  sideLengths.z);
}

void AABB::extend(const AABB& other_box) {
	extend_max(other_box.max);
	extend_min(other_box.min);
}

vec3 AABB::get_center() const {
	return (min + max) * 0.5f;
}

bool AABB::operator==(const AABB& other) {
	return (
		min.x == other.min.x &&
		min.y == other.min.y &&
		min.z == other.min.z &&
		max.x == other.max.x &&
		max.y == other.max.y &&
		max.z == other.max.z
	);
}

bool AABB::operator!=(const AABB& other) {
	return !(*this == other);
}

}