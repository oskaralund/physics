#include "pbd_sandbox.h"
#include "oskgl.h"
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <GLFW/glfw3.h>
#include <iostream>

using namespace glm;

void PBDSandbox::HandleMouseMove(double x, double y)
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

  vertex_manipulator_->MouseMoveCallback(x, y);
}

void PBDSandbox::HandleMouseButton(int button, int action, int mods)
{
  double x, y;
  GetNormalizedMouseCoordinates(&x, &y);
  switch (action)
  {
  case GLFW_PRESS:
    vertex_manipulator_->MouseDownCallback(x, y);
    break;
  case GLFW_RELEASE:
    vertex_manipulator_->MouseUpCallback(x, y);
    break;
  }

}

void PBDSandbox::HandleKey(int key, int scancode, int action, int mods)
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
    case GLFW_KEY_T:
      cloth_->AddForce(0, vec3{10.0f, 0.0f, 0.0f});
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

void PBDSandbox::Reshape(int width, int height)
{
  glViewport(0, 0, width, height);
  float aspect_ratio = static_cast<float>(width)/height;
  camera_.SetAspectRatio(aspect_ratio);
  UpdateProjectionMatrix();
}

PBDSandbox::PBDSandbox()
{
  CreateMatrixBuffer();
  camera_.SetPosition(vec3{0.0f, 0.0f, 1.0f});

  const int N = 35;
  const float side_length = 1.0f;
  /* Construct triangle mesh */
  for (float i = 0; i < N; ++i)
  {
    for (float j = 0; j < N; ++j)
    {
      mesh_.AddVertex(vec3{side_length*i/(N-1), 0.0f, side_length*j/(N-1)});
    }
  }

  for (int i = 0; i < N-1; ++i)
  {
    for (int j = 0; j < N-1; ++j)
    {
      mesh_.AddTriangle(j + i*N, j + i*N+1, j + i*N+N+1);
      mesh_.AddTriangle(j + i*N, j + i*N+N+1, j + i*N+N);
    }
  }

  cloth_ = std::make_unique<PBDFramework>(mesh_);
  vertex_manipulator_ = std::make_unique<VertexManipulator>(cloth_.get(), &camera_);
  vertex_manipulator_->SetPickingForce(10.0f);

  glGenVertexArrays(1, &vao_);
  glBindVertexArray(vao_);
  glGenBuffers(1, &vbo_);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_);
  const std::vector<vec3>& verts = cloth_->GetVertices();
  glBufferData(GL_ARRAY_BUFFER,
               verts.capacity()*sizeof(vec3),
               &verts[0][0],
               GL_DYNAMIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glGenBuffers(1, &element_buffer_);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffer_);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
               mesh_.GetNumTriangles()*3*sizeof(int), &mesh_.GetTriangles()[0][0], GL_DYNAMIC_DRAW);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  mesh_.ComputeVertexNormals();
  glGenBuffers(1, &normal_vbo_);
  glBindBuffer(GL_ARRAY_BUFFER, normal_vbo_);
  const std::vector<vec3>& norms = mesh_.GetVertexNormals();
  glBufferData(GL_ARRAY_BUFFER,
               norms.capacity()*sizeof(vec3),
               &norms[0][0],
               GL_DYNAMIC_DRAW);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

  shader_ = oskgl::CompileShaderProgram("basic_vshader.vsh", "basic_fshader.fsh");
  glUseProgram(shader_);
  glUniformMatrix4fv(0, 1, GL_FALSE, value_ptr(mat4{1.0f}));
}

void PBDSandbox::DrawScene()
{
  UpdateViewMatrix();
  grid_vis_.Draw();
  glUseProgram(shader_);
  glBindBuffer(GL_ARRAY_BUFFER, normal_vbo_);
  mesh_.vertices_ = cloth_->GetVertices();
  mesh_.ComputeVertexNormals();
  const std::vector<vec3>& norms = mesh_.GetVertexNormals();
  glBufferData(GL_ARRAY_BUFFER,
               norms.capacity()*sizeof(vec3),
               &norms[0][0],
               GL_DYNAMIC_DRAW);
  mat4 view = camera_.GetViewMatrix();
  mat3 norm_to_cam = mat3{vec3{view[0]}, vec3{view[1]}, vec3{view[2]}};
  glUniformMatrix4fv(3, 1, GL_FALSE, value_ptr(norm_to_cam));
  vec3 dir_to_light = normalize(vec3{1.0f, 1.0f, 0.0f});
  glUniform3fv(1, 1, value_ptr(dir_to_light));
  vec4 light_intensity = vec4{1.0f, 1.0f, 1.0f, 1.0f};
  glUniform4fv(2, 1, value_ptr(light_intensity));
  glBindBuffer(GL_ARRAY_BUFFER, vbo_);
  const std::vector<vec3>& verts = cloth_->GetVertices();
  glBufferData(GL_ARRAY_BUFFER,
               verts.capacity()*sizeof(vec3),
               &verts[0][0],
               GL_DYNAMIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(vao_);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffer_);
  glDrawElements(GL_TRIANGLES, 3*mesh_.GetNumTriangles(), GL_UNSIGNED_INT, 0);
}

void PBDSandbox::UpdateDynamics(double dt)
{
  MoveCamera(static_cast<float>(dt));
  time_accumulator_ += static_cast<float>(dt);
  while (time_accumulator_ > cloth_->GetTimestep())
  {
    cloth_->Move();
    vertex_manipulator_->Update();
    time_accumulator_ -= cloth_->GetTimestep();
    for (int i = 0; i < cloth_->GetNumVertices(); ++i)
    {
      float y = cloth_->GetVertex(i).y;
      if (y < -0.999f)
      {
        cloth_->Displace(i, vec3{0.0f, -0.999f - y, 0.0f});
        //ridiculous friction
        cloth_->AddForce(i, -cloth_->GetMass(i)*cloth_->GetVelocity(i));
      }
    }
  }
}

void PBDSandbox::CreateMatrixBuffer()
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

void PBDSandbox::UpdateProjectionMatrix()
{
  mat4 projection = camera_.GetProjectionMatrix();
  glBindBuffer(GL_UNIFORM_BUFFER, matrix_buffer_);
  glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(mat4), value_ptr(projection));
  glBindBuffer(GL_UNIFORM_BUFFER, 0);
}

void PBDSandbox::UpdateViewMatrix()
{
  mat4 view = camera_.GetViewMatrix();
  glBindBuffer(GL_UNIFORM_BUFFER, matrix_buffer_);
  glBufferSubData(GL_UNIFORM_BUFFER, sizeof(mat4), sizeof(mat4), value_ptr(view));
  glBindBuffer(GL_UNIFORM_BUFFER, 0);
}

void PBDSandbox::MoveCamera(const float dt)
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
