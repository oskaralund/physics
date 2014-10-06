#include "pbd_framework.h"
#include "oskgl.h"
#include "pbd_kernels.h"
#include <cuda_gl_interop.h>
#include <device_launch_parameters.h>

PBDFramework::PBDFramework()
{
  CreateVertexBuffers();
  InitializeArrays();
}

void PBDFramework::CreateVertexBuffers()
{

  /* Construct grid vertices */
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

  /* Construct index buffer */
  glm::uvec3 grid_elements[grid_size_-1][grid_size_-1][2];
  for (int i = 0; i < grid_size_-1; ++i)
  {
    for (int j = 0; j < grid_size_-1; ++j)
    {
      grid_elements[i][j][0] =
        {i*grid_size_+j, (i+1)*grid_size_+j, (i+1)*grid_size_+j+1};
      grid_elements[i][j][1] =
        {i*grid_size_+j, (i+1)*grid_size_+j+1, i*grid_size_+j+1};
    }
  }

  /* Buffer data */
  glGenBuffers(1, &vertex_vbo_);
  glBindBuffer(GL_ARRAY_BUFFER, vertex_vbo_);
  glBufferData(GL_ARRAY_BUFFER, sizeof(grid_data), grid_data, GL_DYNAMIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glGenBuffers(1, &ibo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(grid_elements), grid_elements, GL_STATIC_READ);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  /* Construct VAOs */
  glGenVertexArrays(1, &vertex_vao_);
  glBindVertexArray(vertex_vao_);
  glBindBuffer(GL_ARRAY_BUFFER, vertex_vbo_);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);

  glBindVertexArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  cudaGraphicsGLRegisterBuffer(&vertices_cuda_,
                               vertex_vbo_,
                               cudaGraphicsMapFlagsWriteDiscard);

  shader_ = oskgl::CompileShaderProgram("cloth.vsh", "cloth.fsh");
}

void PBDFramework::InitializeArrays()
{
  cudaMalloc(&p_, num_vertices_*sizeof(glm::vec3));
  cudaMalloc(&v_, num_vertices_*sizeof(glm::vec3));
  cudaMemset(v_, 0, num_vertices_*sizeof(glm::vec3));
  cudaMalloc(&ext_f_, num_vertices_*sizeof(glm::vec3));
  cudaMemset(ext_f_, 0, num_vertices_*sizeof(glm::vec3));
}

void PBDFramework::Draw()
{
  glUseProgram(shader_);
  glBindVertexArray(vertex_vao_);
  glDrawElements(GL_TRIANGLES, 3*num_triangles_, GL_UNSIGNED_INT, 0);
  glBindVertexArray(0);
  glUseProgram(0);
}

void PBDFramework::Move()
{
  cudaGraphicsMapResources(1, &vertices_cuda_, 0);
  size_t num_bytes;
  cudaGraphicsResourceGetMappedPointer((void**) &x_,
                                       &num_bytes,
                                       vertices_cuda_);
  GuessNewPositions<<<grid_size_, grid_size_>>>(x_,
                                                v_,
                                                ext_f_,
                                                timestep_,
                                                inverse_vertex_mass_,
                                                p_);
  CopyPositions<<<grid_size_, grid_size_>>>(x_, p_);
  SetZero<<<grid_size_, grid_size_>>>(ext_f_);

  cudaGraphicsUnmapResources(1, &vertices_cuda_, 0);
}

PBDFramework::~PBDFramework()
{
  cudaFree(p_);
  cudaFree(v_);
  cudaFree(ext_f_);
}
