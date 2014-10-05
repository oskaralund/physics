#include "pbd_framework.h"
#include <cuda_gl_interop.h>

void PBDFramework::CreateVertexBuffers()
{
  glGenBuffers(2, vertex_vbos_);

  glm::vec3 grid_data[grid_size_][grid_size_];

  for (int i = 0; i < grid_size_; ++i)
  {
    for (int j = 0; j < grid_size_; ++j)
    {
      const float x = static_cast<float>(j)/(grid_size_-1);
      const float z = static_cast<float>(i)/(grid_size_-1);
      grid_data[i][j] = {x, 0.0f, z};
    }
  }

  glBindBuffer(GL_ARRAY_BUFFER, vertex_vbos_[0]);
  glBufferData(GL_ARRAY_BUFFER, sizeof(grid_data), grid_data, GL_DYNAMIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, vertex_vbos_[1]);
  glBufferData(GL_ARRAY_BUFFER, sizeof(grid_data), grid_data, GL_DYNAMIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  cudaGraphicsGLRegisterBuffer(&vertices_cuda_[0],
                               vertex_vbos_[0],
                               cudaGraphicsMapFlagsWriteDiscard);
  cudaGraphicsGLRegisterBuffer(&vertices_cuda_[1],
                               vertex_vbos_[1],
                               cudaGraphicsMapFlagsWriteDiscard);
}
