#include "frame_vis.h"

#include "oskgl.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>

using namespace glm;

FrameVis::FrameVis()
{
  CompileShader();
  ConstructMesh();
}

FrameVis::FrameVis(const bool show_x, const bool show_y, const bool show_z)
  : show_x_axis_(show_x)
  , show_y_axis_(show_y)
  , show_z_axis_(show_z)
{
  CompileShader();
  ConstructMesh();
  num_vertices_ = 2*(show_x_axis_+show_y_axis_+show_z_axis_);
}

void FrameVis::Draw()
{
  glUseProgram(shader_);
  model_ = translate(mat4(1.0f), origo_);
  model_ = model_*orientation_;
  model_ = scale(model_, vec3(scale_));
  glUniformMatrix4fv(0, 1, GL_FALSE, value_ptr(model_));
  glBindVertexArray(vao_);
  glDrawArrays(GL_LINES, 0, num_vertices_);
  glBindVertexArray(0);
}

void FrameVis::ConstructMesh()
{
  assert(sizeof(vec3) == 3*sizeof(float));

  std::vector<vec3> data;
  data.reserve(num_vertices_);

  if (show_x_axis_)
  {
    data.push_back(vec3{0.0f, 0.0f, 0.0f});
    data.push_back(vec3{1.0f, 0.0f, 0.0f});
  }
  if (show_y_axis_)
  {
    data.push_back(vec3{0.0f, 0.0f, 0.0f});
    data.push_back(vec3{0.0f, 1.0f, 0.0f});
  }
  if (show_z_axis_)
  {
    data.push_back(vec3{0.0f, 0.0f, 0.0f});
    data.push_back(vec3{0.0f, 0.0f, 1.0f});
  }

  glGenVertexArrays(1, &vao_);
  glBindVertexArray(vao_);
  GLuint vbo;
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER,
               data.capacity()*sizeof(vec3),
               &data[0],
               GL_STATIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

}

void FrameVis::CompileShader()
{
  shader_ = oskgl::CompileShaderProgram("basic_vshader.vsh", "basic_fshader.fsh");
}

void FrameVis::SetOrigo(const vec3 origo)
{
  origo_ = origo;
}

void FrameVis::SetOrientation(const mat4 orientation)
{
  orientation_ = orientation;
}

void FrameVis::SetScale(const float scale)
{
  scale_ = scale;
}
