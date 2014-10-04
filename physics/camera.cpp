#include "camera.h"

#include <glm/gtc/matrix_transform.hpp>

Camera::Camera()
{
  position_.x        =  0.0f;
  position_.y        =  0.0f;
  position_.z        =  0.0f;
  direction_.x       =  0.0f;
  direction_.y       =  0.0f;
  direction_.z       = -1.0f;
  up_.x              =  0.0f;
  up_.y              =  1.0f;
  up_.z              =  0.0f;
  view_matrix_       =  glm::lookAt(position_, position_+direction_, up_);
}

void Camera::SetPosition(const glm::vec3 &position)
{
  position_    = position;
  view_matrix_ = glm::lookAt(position_, position_+direction_, up_);
}

void Camera::SetDirection(const glm::vec3 &direction)
{
  direction_   = direction;
  view_matrix_ = glm::lookAt(position_, position_+direction_, up_);
}

void Camera::SetUp(const glm::vec3 &up)
{
  up_          = up;
  view_matrix_ = glm::lookAt(position_, position_+direction_, up_);
}
