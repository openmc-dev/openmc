// AABB is short for axis-aligned bounding box. These objects are often used in octrees and BVHs
// Why do we have two bounding box classes? I actually didn't realize openmc had its own class for AABBs until I finished implementing this one
// Regardless, I made this one easier to work with for octrees and BVHs so I'll use this in the partitioner code
// I'll resolve this issue later by trying to merge the two AABB classes into one
#ifndef OPENMC_AABB_H
#define OPENMC_AABB_H

#include "position.h"

namespace openmc {

using vec3 = Position;

struct AABB {
	AABB();
	AABB(const vec3& mi, const vec3& ma);

	void extend_max(const vec3& val);
	void extend_min(const vec3& val);

	void extend(const vec3& pos);
	void extend(const AABB& other_box);

	vec3 get_center() const;
	float surface_area() const;
    bool contains(const vec3& pos) const;

	void reset();

	bool operator==(const AABB& other) const;
	bool operator!=(const AABB& other) const;

    vec3 min;
	vec3 max;
};


}

#endif