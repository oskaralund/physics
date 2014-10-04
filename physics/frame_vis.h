#ifndef FRAME_VIS_H_
#define FRAME_VIS_H_

#include <glload/gl_4_3.h>
#include <glm/glm.hpp>

class FrameVis
{
public:
  FrameVis();
  FrameVis(const bool show_x, const bool show_y, const bool show_z);
  void SetOrigo(const glm::vec3 origo);
  void SetOrientation(const glm::mat4 orientation);
  void SetScale(const float scale);
  void Draw();

private:
  void ConstructMesh();
  void CompileShader();

  bool show_x_axis_ = true;
  bool show_y_axis_ = true;
  bool show_z_axis_ = true;
  int num_vertices_ = 0;
  float scale_ = 1.0f;
  glm::vec3 origo_{0.0f,0.0f,0.0f};
  glm::mat4 orientation_{1.0f};
  GLuint shader_;
  GLuint vao_;
  glm::mat4 model_{1.0f};
};

#endif
