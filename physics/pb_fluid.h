#ifndef PB_FLUID_H_
#define PB_FLUID_H_

#include <glm/glm.hpp>

class PBFluid
{
public:
  PBFluid();
  void Move();
private:
  void Construct3DParticleGrid();

  timestep_ = 1e-2f;
  std::vector<glm::vec3> vertices_;
  std::vector<glm::vec3> new_vertices_;
  std::vector<glm::vec3> velocities_;
  std::vector<glm::vec3> external_forces_;
};

#endif
