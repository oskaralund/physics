#ifndef ROD_INTERFACE_H_
#define ROD_INTERFACE_H_

#include <glm/glm.hpp>

class RodInterface
{
public:
  virtual void Move() = 0;
  virtual void SetTimestep(float) = 0;
  virtual void Displace(int, const glm::vec3&) = 0;
  virtual int GetNumVertices() const = 0;
  virtual int GetNumEdges() const = 0;
  virtual float GetTimestep() const = 0;
  virtual float GetEdgeRadius(int) const = 0;
  virtual float GetEdgeLength(int) const = 0;
  virtual glm::vec3 GetVertex(int) const = 0;
  virtual glm::mat4 GetMaterialFrame(int) const = 0;
};

#endif
