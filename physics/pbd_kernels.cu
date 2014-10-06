#include <glm/glm.hpp>

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
