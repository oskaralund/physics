#ifndef CYLINDER_VIS_H_
#define CYLINDER_VIS_H_

#include <memory>

#include <glload/gl_4_3.h>
#include <glm/glm.hpp>

class CylinderVis
{
public:

  CylinderVis(const int num_vertices);
  void Draw();
  void SetRadius(const float r);
  void SetStartPoint(const glm::vec3 start);
  void SetEndPoint(const glm::vec3 end);
  void SetOrientation(const glm::mat4 orientation);
  void SetLength(const float length);

private:
  void ConstructMesh();
  void ComputeModelMatrix();
  void CompileShader();

  int num_vertices_;
  GLuint shader_;
  GLuint vao_;
  float radius_ = 1.0f;
  float length_ = 1.0f;
  glm::vec3 start_{0.0f,0.0f,0.0f};
  glm::vec3 end_{0.0f,0.0f,1.0f};
  glm::mat4 model_{1.0f};
  glm::mat4 orientation_{1.0f};

};

#endif
