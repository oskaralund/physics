#ifndef FPS_CAMERA_H_
#define FPS_CAMERA_H_

#include "camera_interface.h"
#include <glm/glm.hpp>

class FPSCamera : public CameraInterface
{
public:
  void Yaw(const float angle);
  void Pitch(const float angle);
  inline void FPSCamera::MoveForward(const float distance)
  {
    Move(distance*direction_);
  }
  inline void FPSCamera::MoveBack(const float distance)
  {
    Move(-distance*direction_);
  }
  inline void FPSCamera::MoveLeft(const float distance)
  {
    Move(-distance*right_);
  }
  inline void FPSCamera::MoveRight(const float distance)
  {
    Move(distance*right_);
  }
  inline void FPSCamera::MoveUp(const float distance)
  {
    Move(distance*up_);
  }
  inline void FPSCamera::MoveDown(const float distance)
  {
    Move(-distance*up_);
  }
  void SetPosition(const glm::vec3 pos) { position_ = pos; }
  void SetFieldOfView(const float f) { fov_ = f; }
  void SetZNear(const float z_near) { z_near_ = z_near; }
  void SetZFar(const float z_far) { z_far_ = z_far; }
  void SetAspectRatio(const float aspect) { aspect_ratio_ = aspect; }
  glm::vec3 GetPosition() const { return position_; }
  glm::mat4 GetViewMatrix() const override;
  glm::mat4 GetProjectionMatrix() const override;

private:
  inline void Move(const glm::vec3 d) { position_ += d; }
  float pitch_ = 0.0f;
  float fov_ = 90.0f;
  float z_near_ = 0.01f;
  float z_far_ = 100.0f;
  float aspect_ratio_ = 16.0f/9.0f;
  glm::vec3 position_{0.0f};
  glm::vec3 direction_{0, 0, -1};
  glm::vec3 right_{1, 0, 0};
  glm::vec3 up_{0, 1, 0};
};

#endif
