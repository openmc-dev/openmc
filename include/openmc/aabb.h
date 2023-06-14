// AABB is short for axis-aligned bounding box. These objects are often used in octrees and BVHs
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

	void reset();

	bool operator==(const AABB& other);
	bool operator!=(const AABB& other);

    vec3 min;
	vec3 max;
};


}

#endif