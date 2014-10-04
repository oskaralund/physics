#ifndef DISTANCE_SEGMENT_SEGMENT_H
#define DISTANCE_SEGMENT_SEGMENT_H

#include <glm/glm.hpp>

inline float DistanceLinesegLineseg(const glm::vec3 p0, const glm::vec3 p1, const glm::vec3 q0, const glm::vec3 q1, float* closest_point_seg0, float* closest_point_seg1);

#include "distance_lineseg_lineseg.inl"

#endif
