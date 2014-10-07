#ifndef DISTANCE_POINT_RAY_H_
#define DISTANCE_POINT_RAY_H_

#include <glm/glm.hpp>

/* Computes the distance of a point p to the closest point on a ray.
 * ray_direction must be normalized.
 */
inline float DistancePointRay(const glm::vec3& p,
                              const glm::vec3& ray_origin,
                              const glm::vec3& ray_direction,
                              glm::vec3* closest_point);

#include "distance_point_ray.inl"

#endif
