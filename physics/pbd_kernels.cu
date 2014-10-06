#include "pbd_kernels.h"

__device__ int grid_size = 10;

__device__ void ProjectLengthConstraint(glm::vec3* p, int i, int j, float d)
{
  glm::vec3 p1 = p[i];
  glm::vec3 p2 = p[j];
  p[i] -= 0.5f*(glm::length(p1 - p2) - d)*glm::normalize(p1-p2);
  p[j] += 0.5f*(glm::length(p1 - p2) - d)*glm::normalize(p1-p2);
}

__global__ void GuessNewPositions(glm::vec3* x,
                                  glm::vec3* v,
                                  glm::vec3* f,
                                  const float dt,
                                  const float w,
                                  glm::vec3* p)
{
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  v[i] += dt*w*f[i];
  p[i] = x[i] + dt*v[i];
}

__global__ void UpdateVelocities(glm::vec3* x,
                                 glm::vec3* p,
                                 glm::vec3* v,
                                 const float dt)
{
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  v[i] += (p[i] - x[i])/dt;
}

__global__ void ProjectLengthConstraints(glm::vec3* p)
{
  const float side_length = 1.0f/(grid_size-1);
  const float diag_length = glm::sqrt(2*side_length*side_length);
  if (blockDim.x & 1)
  {
    int i = (gridDim.x+1)*threadIdx.x + blockIdx.x;
    int j = i+1;
    ProjectLengthConstraint(p, i, j, side_length);
  }
  __syncthreads();
  if (!(blockDim.x & 1))
  {
    int i = (gridDim.x+1)*threadIdx.x + blockIdx.x;
    int j = i+1;
    ProjectLengthConstraint(p, i, j, side_length);
  }
  __syncthreads();
  if (blockDim.x & 1)
  {
    int i = (gridDim.x+1)*blockIdx.x + threadIdx.x;
    int j = i+gridDim.x+1;
    ProjectLengthConstraint(p, i, j, side_length);
  }
  __syncthreads();
  if (!(blockDim.x & 1))
  {
    int i = (gridDim.x+1)*blockIdx.x + threadIdx.x;
    int j = i+gridDim.x+1;
    ProjectLengthConstraint(p, i, j, side_length);
  }
  //__syncthreads();
  //if (threadIdx.x < grid_size-1)
  //{
  //  int i = (gridDim.x+1)*threadIdx.x + blockIdx.x;
  //  int j = i+gridDim.x+2;
  //  ProjectLengthConstraint(p, i, j, diag_length);
  //}
  //__syncthreads();
  //if (threadIdx.x < grid_size-1)
  //{
  //  int i = (gridDim.x+1)*threadIdx.x + blockIdx.x + 1;
  //  int j = i+gridDim.x;
  //  ProjectLengthConstraint(p, i, j, diag_length);
  //}
}

__global__ void SetElement(glm::vec3* array, const int i, const glm::vec3 value)
{
  array[i] = value;
}

__global__ void SetZero(glm::vec3* array)
{
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  array[i] = {0.0f, 0.0f, 0.0f};
}

__global__ void CopyPositions(glm::vec3* x,
                              glm::vec3* p)
{
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  x[i] = p[i];
}
