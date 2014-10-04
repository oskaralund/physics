#ifndef DYNAMICAL_SYSTEM_INTERFACE_H_
#define DYNAMICAL_SYSTEM_INTERFACE_H_

#include <vector>
#include <glm/glm.hpp>

class DynamicalSystemInterface
{
public:
  virtual void Displace(int, const glm::vec3&) = 0;
  virtual void AddForce(int, const glm::vec3&) = 0;
  virtual float GetVertexRadius(int) const = 0;
  virtual glm::vec3 GetVertex(int) const = 0;
  virtual glm::vec3 GetVelocity(int) const = 0;
  virtual const std::vector<glm::vec3>& GetVertices() const = 0;
};

#endif
