#include "fps_camera.h"

#include <glm/gtc/constants.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <GLFW/glfw3.h>

using namespace glm;

mat4 FPSCamera::GetViewMatrix() const
{
  return lookAt(position_,
                position_ + direction_,
                up_);
}

mat4 FPSCamera::GetProjectionMatrix() const
{
  return perspective(pi<float>()/2, aspect_ratio_, z_near_, z_far_);
}

void FPSCamera::Yaw(const float angle)
{
  direction_ = rotate(direction_, angle, -up_);
  right_ = rotate(right_, angle, -up_);
}

void FPSCamera::Pitch(const float angle)
{
  float delta_angle = angle;
  if (angle < 0 && pitch_ < -0.95f*half_pi<float>())
    delta_angle = 0.0f;
  if (angle > 0 && pitch_ > 0.95f*half_pi<float>())
    delta_angle = 0.0f;

  pitch_ += delta_angle;
  direction_ = rotate(direction_, delta_angle, right_);
}
