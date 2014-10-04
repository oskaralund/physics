#ifndef CAMERA_INTERFACE_H_
#define CAMERA_INTERFACE_H_

#include <glm/glm.hpp>

class CameraInterface
{
public:
  virtual glm::mat4 GetViewMatrix() const = 0;
  virtual glm::mat4 GetProjectionMatrix() const = 0;
};

#endif
