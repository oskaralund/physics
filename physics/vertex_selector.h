#ifndef VERTEX_SELECTOR_H_
#define VERTEX_SELECTOR_H_

#include "vertex_selector_interface.h"

class VertexSelector : public VertexSelectorInterface
{
public:
  VertexSelector(DynamicalSystemInterface* system) : dynamical_system_(system) { }
  void SetDynamicalSystem(DynamicalSystemInterface*) override;
  int GetVertexClosestToRay(const glm::vec3& ray_origin,
                            const glm::vec3& ray_direction) const override;
private:
  DynamicalSystemInterface* dynamical_system_ = nullptr;
};

#endif
