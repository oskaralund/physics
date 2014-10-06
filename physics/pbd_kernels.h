/*
The kernels correspond to steps in the PBD algorithm
as described in Müller et al. (2006). Parameter names
have been chosen to reflect the names in the algorithm:

N = number of vertices.
x_i = vertex i.
v_i = velocity i.
f_ext_i = external force acting on vertex i.
w_i = inverse mass of vertex i.

M = number of constraints (excluding collision constraints).
M_coll = number of collision constraints.
C_i = constraint i.

loop
  forall vertices i do v_i <- v_i + dt*w_i*f_ext_i
  dampVelocities(v_1, ..., v_N)
  forall vertices i do p_i <- x_i + dt*v_i
  forall vertices i do generateCollisionConstraints(x_i -> p_i)
  loop solverIterations times
    projectConstraints(C_1, ..., C_{M+M_coll}, p_1, ..., p_N)
  endloop
  forall vertices i
    v_i <- (p_i - x_i)/dt
    x_i <- p_i
  endfor
endloop
*/

__global__ void GuessNewPositions(glm::vec3* x,   //Input: Positions
                                  glm::vec3* v,   //Input: Velocities
                                  glm::vec3* f,   //Input: External forces
                                  const float dt, //Input: Timestep
                                  const float w,  //Input: Inverse vertex mass
                                  glm::vec3* p);  //Output: New positions

__global__ void ProjectLengthConstraints(glm::vec3* p);

__global__ void UpdateVelocities(glm::vec3* x,    //Input: Old positions
                                 glm::vec3* p,    //Input: New positions
                                 glm::vec3* v,    //Input/Output: Velocities
                                 const float dt); //Input: Timestep

/* Set element i of array to value */
__global__ void SetElement(glm::vec3* array, int i, glm::vec3 value);

/* Set all elements in array to zero */
__global__ void SetZero(glm::vec3* array);

/* Copy p into x */
__global__ void CopyPositions(glm::vec3* x,
                              glm::vec3* p);

/* Convert matrix indices to a vector index (assuming row major order) */
__device__ int ToVectorIndex(int i, int j, int ld)
{
  return i*ld + j;
}
