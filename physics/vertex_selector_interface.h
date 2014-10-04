#ifndef VERTEX_SELECTOR_INTERFACE_H_
#define VERTEX_SELECTOR_INTERFACE_H_

#include <glm/glm.hpp>
#include "dynamical_system_interface.h"

class VertexSelectorInterface
{
public:
  virtual void SetDynamicalSystem(DynamicalSystemInterface*) = 0;
  virtual int GetVertexClosestToRay(const glm::vec3& ray_origin,
                                    const glm::vec3& ray_direction) const = 0;
};

#endif
