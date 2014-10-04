#ifndef GRID_VIS_H_
#define GRID_VIS_H_

#include <glload/gl_4_3.h>
#include <glm/glm.hpp>

class GridVis
{
public:
  GridVis(const int grid_size,
          const float square_size,
          const glm::vec3 origin,
          const glm::vec3 normal);
  void Draw();
  void SetGridSize(const int);
  void SetSquareSize(const float);
  void SetOrigin(const glm::vec3);
  void SetNormal(const glm::vec3);

private:
  void ConstructMesh();
  void ComputeModelMatrix();
  void CompileShader();

  GLuint shader_;
  GLuint vao_;
  int grid_size_ = 0;
  float square_size_ = 0;
  glm::vec3 origin_{0.0f, 0.0f, 0.0f};
  glm::vec3 normal_{0.0f, 1.0f, 0.0f};
  glm::mat4 model_{1.0f};
};

#endif
