#ifndef BASIC_ROD_H_
#define BASIC_ROD_H_

#include <vector>

#include "rod_interface.h"
#include "dynamical_system_interface.h"

class BasicRod : public RodInterface, public DynamicalSystemInterface
{
public:
  BasicRod(const float rod_length,
           const int num_edges,
           const float radius,
           const float density);

  void  Move() override;
  void  AddForce(int i, const glm::vec3& f) override { external_forces_[i] += f; }
  void  SetTimestep(float h) override { timestep_ = h; }
  void  Displace(int i, const glm::vec3& d) override;
  int   GetNumVertices() const override { return kNumVertices; }
  int   GetNumEdges() const override { return kNumEdges; }
  float GetTimestep() const override { return timestep_; }
  float GetVertexRadius(int i) const override;
  float GetEdgeRadius(int i) const override;
  float GetEdgeLength(int i) const override { return glm::length(edges_[i]); }
  glm::vec3 GetVertex(int i) const override { return vertices_[i]; }
  glm::vec3 GetVelocity(int i) const override { return velocities_[i]; }
  glm::mat4 GetMaterialFrame(int i) const override;
  const std::vector<glm::vec3>& GetVertices() const override { return vertices_; }

private:
  enum RelativePosition { PREV, CURR, NEXT };
  void Init();
  void IntegrateForces();
  void EnforceLengthConstraints();
  void ComputeEdges();
  void ComputeBinorms();
  void ComputeConstraints();
  void ComputeConstraintMatrix();
  void ComputeReferenceTwist();
  void ComputeMaterialFrames();
  void ComputeInternalForces();
  void UpdateVelocities();
  void UpdateReferenceFrames();
  void ResetExternalForces();

  float GetTwist(int i) const;
  float GetStretch(int i) const;
  float GetTwistEnergyAngleDerivative(int i) const;
  float GetBendEnergyAngleDerivative(int i) const;
  float GetStretchStiffness(int i) const;
  float GetBendStiffness(int i) const;
  float GetTwistStiffness(int i) const;
  float GetEdgeCrossSectionArea(int i) const;
  float GetVertexCrossSectionArea(int i) const;
  float GetCrossSectionalInertia(int i) const;

  bool viscous_ = true;

  glm::vec2 GetCurvature(int i) const;
  glm::vec2 GetRestCurvature(int i) const;
  glm::vec2 GetCurvatureAngleJacobian(int i, int j) const;
  glm::vec2 GetCurvatureAngleJacobianJacobian(int i, int j) const;
  glm::vec3 GetTwistVertexGradient(int i, int j) const;
  glm::vec3 GetTwistEdgeGradient(int i, int j) const;
  glm::vec3 GetBendEnergyVertexGradient(int i) const;
  glm::vec3 GetTwistEnergyVertexGradient(int i) const;
  glm::vec3 GetStretchEnergyVertexGradient(int i) const;
  glm::vec3 GetStretchVertexGradient(int i, int j) const;
  glm::vec3 GetStretchEdgeGradient(int i) const;

  glm::mat2x3 GetCurvatureVertexJacobian(int i, int j) const;
  glm::mat2x3 GetCurvatureEdgeJacobian(int i, int j) const;

  const int kNumEdges;
  const int kNumVertices = kNumEdges+1;

  const float kRodLength;
  const float kEdgeLength = kRodLength/kNumEdges;
  const float kEdgeLengthSquared = kEdgeLength*kEdgeLength;
  const float kEdgeLengthInv = 1/kEdgeLength;
  const float kEdgeLengthInvSquared = 1/(kEdgeLength*kEdgeLength);
  const float kRadius;
  const float kDensity;
  const float kVertexMass;
  const float kVertexMassInv = 1/kVertexMass;
  const float kAngleTolerance = 2.0f*3.1415f/1000.0f;
  float timestep_ = 1e-3f;
  float young_modulus_ = 1e3f;
  float shear_modulus_ = 1e3f;
  float damping_ = 0.005f;
  float viscosity_ = 5.5f;

  glm::vec3 gravity_{0.0f, -9.82f, 0.0f};

  std::vector<float> constraints_;
  std::vector<float> constr_diagonal_;
  std::vector<float> constr_off_diagonal_;
  std::vector<float> energy_angle_gradient_;
  std::vector<float> energy_angle_hessian_diag_;
  std::vector<float> energy_angle_hessian_off_diag_;
  std::vector<float> reference_twist_;
  std::vector<float> material_frame_angles_;
  std::vector<float> angle_velocities_;
  std::vector<float> rest_stretch_;
  std::vector<glm::vec2> rest_curvatures_;
  std::vector<glm::vec3> vertices_;
  std::vector<glm::vec3> prev_vertices_;
  std::vector<glm::vec3> velocities_;
  std::vector<glm::vec3> edges_;
  std::vector<glm::vec3> prev_edges_;
  std::vector<glm::vec3> binorms_;
  std::vector<glm::vec3> external_forces_;
  std::vector<glm::vec3> internal_forces_;
  std::vector<std::pair<glm::vec3, glm::vec3>> reference_frames_;
  std::vector<std::pair<glm::vec3, glm::vec3>> material_frames_;
};

#endif
