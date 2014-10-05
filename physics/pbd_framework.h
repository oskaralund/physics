#ifndef PBD_FRAMEWORK_H_
#define PBD_FRAMEWORK_H_

#include <cuda_runtime.h>
#include <glload/gl_4_3.h>
#include <glload/gl_load.h>
#include <glm/glm.hpp>

class PBDFramework
{
public:

private:
  void CreateVertexBuffers();
  static const int grid_size_ = 10;
  GLuint vertex_vbos_[2];
  glm::vec3* vertices_ = nullptr;
  struct cudaGraphicsResource* vertices_cuda_[2];
};

#endif
