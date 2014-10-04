inline float DistancePointRay(const glm::vec3& p,
                              const glm::vec3& ray_origin,
                              const glm::vec3& ray_direction,
                              glm::vec3* closest_point)
{
  float proj = glm::dot(p - ray_origin, ray_direction);
  if (proj > 0.0f)
    *closest_point = ray_origin + proj*ray_direction;
  else
    *closest_point = ray_origin;

  return glm::length(*closest_point - p);
}
