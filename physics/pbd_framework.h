#ifndef PBD_FRAMEWORK_H_
#define PBD_FRAMEWORK_H_

#include <cuda_runtime.h>
#include <glload/gl_4_3.h>
#include <glload/gl_load.h>
#include <glm/glm.hpp>

class PBDFramework
{
public:
  PBDFramework();
  ~PBDFramework();
  void Draw();
  void Move();
  float GetTimestep() const { return timestep_; }
private:
  void CreateVertexBuffers();
  void InitializeArrays();

  static const int grid_size_ = 4;
  static const int num_triangles_ = 2*(grid_size_-1)*(grid_size_-1);
  static const int num_vertices_ = grid_size_*grid_size_;
  float mass_ = 0.5f;
  float timestep_ = 1e-2f;
  float inverse_vertex_mass_ = num_vertices_/mass_;
  glm::vec3* x_ = nullptr; //Positions
  glm::vec3* p_ = nullptr; //New guessed positions
  glm::vec3* v_ = nullptr; //Velocities
  glm::vec3* ext_f_ = nullptr; //External forces
  struct cudaGraphicsResource* vertices_cuda_;
  GLuint vertex_vbo_;
  GLuint vertex_vao_;
  GLuint ibo;
  GLuint shader_;
};

#endif
