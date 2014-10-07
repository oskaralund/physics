#include "cylinder_vis.h"
#include "oskgl.h"

#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <glm/gtx/quaternion.hpp>

#include <vector>

using namespace glm;

CylinderVis::CylinderVis(const int num_vertices) : num_vertices_(num_vertices)
{
  CompileShader();
  ConstructMesh();
}

void CylinderVis::Draw()
{
  ComputeModelMatrix();
  glUseProgram(shader_);
  glUniformMatrix4fv(0, 1, GL_FALSE, value_ptr(model_));
  glBindVertexArray(vao_);
  glDrawArrays(GL_TRIANGLE_STRIP, 0, num_vertices_);
  glBindVertexArray(0);
}

void CylinderVis::CompileShader()
{
  shader_ = oskgl::CompileShaderProgram("shaders/cylinder.vsh",
                                        "shaders/cylinder.fsh");
}

void CylinderVis::ComputeModelMatrix()
{
  model_ = translate(mat4(1.0f), start_);
  model_ = model_*orientation_;
  model_ = scale(model_, vec3{radius_, radius_, length_});
}

void CylinderVis::ConstructMesh()
{
  assert(sizeof(vec3) == 3*sizeof(float));

  /* Construct vertex data */
  std::vector<vec3> data;
  data.reserve(num_vertices_*2);
  for (int i = 0; i < num_vertices_; i = i+2)
  {
    float theta;
    theta = 2*glm::pi<float>()*i/(num_vertices_-2);
    data.push_back(vec3{cosf(theta), sinf(theta), 0.0f});
    data.push_back(vec3{0.0f, 0.0f, glm::abs(sinf(2.5f*theta))});

    theta = 2*glm::pi<float>()*(i+0.5f)/(num_vertices_-2);
    data.push_back(vec3{cosf(theta), sinf(theta), 1.0f});
    data.push_back(vec3{0.0f, 0.0f, glm::abs(sinf(2.5f*theta))});
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
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 2*sizeof(vec3), 0);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 2*sizeof(vec3), (GLvoid*)sizeof(vec3));
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

void CylinderVis::SetRadius(const float r)
{
  radius_ = r;
}

void CylinderVis::SetStartPoint(const vec3 start)
{
  start_ = start;
}

void CylinderVis::SetEndPoint(const vec3 end)
{
  end_ = end;
}

void CylinderVis::SetOrientation(const mat4 orientation)
{
  orientation_ = orientation;
}

void CylinderVis::SetLength(const float length)
{
  length_ = length;
}
