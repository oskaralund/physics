#ifndef PBD_FRAMEWORK_H_
#define PBD_FRAMEWORK_H_

#include "dynamical_system_interface.h"
#include "triangle_mesh.h"
#include <memory>
#include <array>

class PBDFramework : public DynamicalSystemInterface
{
public:
  PBDFramework();
  PBDFramework(const TriangleMesh& mesh);
  void Move();
  void Displace(int, const glm::vec3&) override;
  void AddForce(int, const glm::vec3&) override;
  float GetVertexRadius(int) const override { return vertex_radius_; }
  float GetMass(int i) const { return 1.0f/inverse_masses_[i]; }
  glm::vec3 GetVertex(int i) const override { return vertices_[i]; }
  glm::vec3 GetVelocity(int i) const override { return velocities_[i]; }
  const std::vector<glm::vec3>& GetVertices() const override { return vertices_; }
  size_t GetNumVertices() const { return num_vertices_; }
  float GetTimestep() const { return timestep_; }

  void AddLengthConstraint(int i, int j, float l);
  void AddBendConstraint(std::array<int, 4> dep_vars,
                         const float b);
  void AddBendConstraint(std::pair<glm::uvec3, glm::uvec3>,
                         const float angle);
  void LoadTriangleMesh(const TriangleMesh& mesh);
  void Resize(const int i);

private:
  void UpdateVelocities();
  void ResetExternalForces();
  void ProjectConstraints();
  void ProjectLengthConstraints();
  void ProjectBendConstraints();
  void IntegrateExternalForces();
  void DampVelocities();
  glm::vec3 GetCenterOfMass();
  glm::vec3 GetCenterOfMassVelocity();
  glm::vec3 GetAngularMomentum();

  /* LENGTH CONSTRAINT */
  int num_length_constraints_;
  float c_length_stiffness_ = 1.0f;
  std::vector<std::pair<int, int>> c_length_dep_vars_;
  std::vector<float> c_length_lengths_;
  /* END LENGTH CONSTRAINT */

  /* BEND CONSTRAINT */
  int num_bend_constraints_;
  const float c_bend_stiffness_ = 0.9f;
  std::vector<float> c_bend_angles_;
  std::vector<std::array<int, 4>> c_bend_dep_vars_;
  /* END BEND CONSTRAINT */

  size_t num_vertices_;
  int num_solver_iterations_ = 8;
  float timestep_ = 1.0f/100.0f;
  float vertex_radius_ = 0.05f;
  float damping_ = 0.01f;
  float mass_ = 0.5f;
  glm::vec3 gravity_{0.0f, -9.82f, 0.0f};
  std::vector<float> inverse_masses_;
  std::vector<glm::vec3> vertices_;
  std::vector<glm::vec3> new_vertices_;
  std::vector<glm::vec3> velocities_;
  std::vector<glm::vec3> external_forces_;
};

#endif
