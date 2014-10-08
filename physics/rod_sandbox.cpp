#include "rod_sandbox.h"
#include "oskgl.h"
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <GLFW/glfw3.h>
#include <iostream>

using namespace glm;

RodSandbox::RodSandbox()
{
  CreateMatrixBuffer();
  camera_.SetPosition(vec3{0.0f, 0.0f, 0.2f});
  rod_vis_.DrawMaterialFrames(true);
  vertex_manipulator_.SetPickingForce(5.0f);
}

void RodSandbox::UpdateDynamics(double dt)
{
  MoveCamera(static_cast<float>(dt));
  time_accumulator_ += static_cast<float>(dt);
  while (time_accumulator_ > rod_.GetTimestep())
  {
    rod_.Move();
    for (int i = 0; i < rod_.GetNumVertices(); ++i)
    {
      if (rod_.GetVertex(i)[1] < floor_level_+rod_.GetEdgeRadius(min(i, rod_.GetNumEdges()-1)))
        rod_.Displace(i, vec3{0.0f, -(rod_.GetVertex(i)[1]-floor_level_-rod_.GetEdgeRadius(min(i, rod_.GetNumEdges()-1))), 0.0f});
    }

    vertex_manipulator_.Update();
    time_accumulator_ -= rod_.GetTimestep();
  }
}

void RodSandbox::DrawScene()
{
  UpdateViewMatrix();
  rod_vis_.Draw();
  grid_vis_.Draw();
}

void RodSandbox::HandleMouseMove(double x, double y)
{
  AbsoluteToNormalizedCoordinates(&x, &y);
  float dx = static_cast<float>(x - prev_x_);
  float dy = static_cast<float>(y - prev_y_);
  const float kPi = pi<float>();
  float yaw_angle = camera_sensitivity_*kPi*dx;
  float pitch_angle = camera_sensitivity_*kPi*dy;
  camera_.Yaw(yaw_angle);
  camera_.Pitch(pitch_angle);

  CenterCursor();
  prev_x_ = 0.0f;
  prev_y_ = 0.0f;

  vertex_manipulator_
    .MouseMoveCallback(x, y);
}

void RodSandbox::HandleMouseButton(int button, int action, int mods)
{
  double x, y;
  GetNormalizedMouseCoordinates(&x, &y);
  switch (action)
  {
  case GLFW_PRESS:
    vertex_manipulator_.MouseDownCallback(x, y);
    break;
  case GLFW_RELEASE:
    vertex_manipulator_.MouseUpCallback(x, y);
    break;
  }

}

void RodSandbox::HandleKey(int key, int scancode, int action, int mods)
{
  switch (action)
  {
  case GLFW_PRESS:
    switch (key)
    {
    case GLFW_KEY_W:
      move_forward_ = true;
      break;
    case GLFW_KEY_A:
      move_left_ = true;
      break;
    case GLFW_KEY_S:
      move_back_ = true;
      break;
    case GLFW_KEY_D:
      move_right_ = true;
      break;
    case GLFW_KEY_SPACE:
      move_up_ = true;
      break;
    case GLFW_KEY_LEFT_CONTROL:
      move_down_ = true;
      break;
    case GLFW_KEY_ESCAPE:
      glfwSetWindowShouldClose(GetGLFWWindowPointer(), 1);
      break;
    }
    break;

  case GLFW_RELEASE:
    switch (key)
    {
    case GLFW_KEY_W:
      move_forward_ = false;
      break;
    case GLFW_KEY_A:
      move_left_ = false;
      break;
    case GLFW_KEY_S:
      move_back_ = false;
      break;
    case GLFW_KEY_D:
      move_right_ = false;
      break;
    case GLFW_KEY_SPACE:
      move_up_ = false;
      break;
    case GLFW_KEY_LEFT_CONTROL:
      move_down_ = false;
      break;
    }
    break;
  }
}

void RodSandbox::Reshape(int width, int height)
{
  glViewport(0, 0, width, height);
  float aspect_ratio = static_cast<float>(width)/height;
  camera_.SetAspectRatio(aspect_ratio);
  UpdateProjectionMatrix();
}

void RodSandbox::CreateMatrixBuffer()
{
  mat4 projection = camera_.GetProjectionMatrix();
  mat4 view = camera_.GetViewMatrix();
  glGenBuffers(1, &matrix_buffer_);
  glBindBuffer(GL_UNIFORM_BUFFER, matrix_buffer_);
  glBufferData(GL_UNIFORM_BUFFER, 2*sizeof(mat4), NULL, GL_STREAM_DRAW);
  glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(mat4), value_ptr(projection));
  glBufferSubData(GL_UNIFORM_BUFFER, sizeof(mat4), sizeof(mat4), value_ptr(view));
  glBindBuffer(GL_UNIFORM_BUFFER, 0);
  glBindBufferBase(GL_UNIFORM_BUFFER, 0, matrix_buffer_);
}

void RodSandbox::UpdateProjectionMatrix()
{
  mat4 projection = camera_.GetProjectionMatrix();
  glBindBuffer(GL_UNIFORM_BUFFER, matrix_buffer_);
  glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(mat4), value_ptr(projection));
  glBindBuffer(GL_UNIFORM_BUFFER, 0);
}

void RodSandbox::UpdateViewMatrix()
{
  mat4 view = camera_.GetViewMatrix();
  glBindBuffer(GL_UNIFORM_BUFFER, matrix_buffer_);
  glBufferSubData(GL_UNIFORM_BUFFER, sizeof(mat4), sizeof(mat4), value_ptr(view));
  glBindBuffer(GL_UNIFORM_BUFFER, 0);
}

void RodSandbox::MoveCamera(const float dt)
{
  if (move_forward_)
    camera_.MoveForward(dt*camera_sensitivity_);
  if (move_back_)
    camera_.MoveBack(dt*camera_sensitivity_);
  if (move_left_)
    camera_.MoveLeft(dt*camera_sensitivity_);
  if (move_right_)
    camera_.MoveRight(dt*camera_sensitivity_);
  if (move_up_)
    camera_.MoveUp(dt*camera_sensitivity_);
  if (move_down_)
    camera_.MoveDown(dt*camera_sensitivity_);
}
