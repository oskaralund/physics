#include "vertex_selector.h"

#include <glm/gtx/intersect.hpp>

using namespace glm;

void VertexSelector::SetDynamicalSystem(DynamicalSystemInterface* system)
{
  dynamical_system_ = system;
}

/* Returns the index of the vertex closest to a ray */
int VertexSelector::GetVertexClosestToRay(const vec3& ray_origin,
                                          const vec3& ray_direction) const
{
  auto vertices = dynamical_system_->GetVertices();
  float min_distance = FLT_MAX;
  float distance;
  int index = -1;
  for (int i = 0; i < vertices.size(); ++i)
  {
    float radius_squared =
      dynamical_system_->GetVertexRadius(i)*dynamical_system_->GetVertexRadius(i);
    bool hit = intersectRaySphere(ray_origin,
                                  ray_direction,
                                  vertices[i],
                                  radius_squared,
                                  distance);
    if (hit && distance < min_distance)
    {
      min_distance = distance;
      index = i;
    }
  }

  return index;
}
