#include "grid_vis.h"
#include "oskgl.h"
#include <glm/gtx/vector_angle.hpp>
#include <glm/gtx/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>

using namespace glm;

GridVis::GridVis(const int grid_size,
                 const float square_size,
                 const vec3 origin,
                 const vec3 normal)
  : grid_size_(grid_size)
  , square_size_(square_size)
  , origin_(origin)
  , normal_(normal)
{
  ConstructMesh();
  CompileShader();

  float half_side = grid_size*square_size/2;
  model_ = translate(mat4{1.0f}, origin);
  float theta = angle(normal_, vec3{0.0f, 1.0f, 0.0f});
  if (theta > 1e-2f*pi<float>())
    model_ = rotate(model_, theta, cross(normal_, vec3{0.0f, 1.0f, 0.0f}));

  model_ = scale(model_, vec3{square_size, 1.0f, square_size});

  glUseProgram(shader_);
  glUniformMatrix4fv(0, 1, GL_FALSE, value_ptr(model_));
  glUniform1i(1, grid_size_);
  glUniform1f(2, square_size_);
  glUseProgram(0);
}

void GridVis::CompileShader()
{
  shader_ = oskgl::CompileShaderProgram("grid_vis.vsh", "grid_vis.fsh");
}

void GridVis::Draw()
{
  glUseProgram(shader_);
  glBindVertexArray(vao_);
  glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, grid_size_*grid_size_);
  glBindVertexArray(0);
  glUseProgram(0);
}

void GridVis::ConstructMesh()
{
  assert(sizeof(vec3) == 3*sizeof(float));

  /* Construct vertex data */
  std::vector<vec3> data(4);
  data[0] = vec3{0.0f, 0.0f, 0.0f};
  data[1] = vec3{1.0f, 0.0f, 0.0f};
  data[2] = vec3{0.0f, 0.0f, 1.0f};
  data[3] = vec3{1.0f, 0.0f, 1.0f};

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
