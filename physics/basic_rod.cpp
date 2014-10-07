#include "basic_rod.h"

#include <algorithm>
#include <iostream>

#include <glm/gtc/constants.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <glm/gtx/norm.hpp>
#include "cblas.h"
#define HAVE_LAPACK_CONFIG_H
#define LAPACK_COMPLEX_STRUCTURE
#include "lapacke.h"

using namespace glm;

namespace {
  /*
  This function computes the smallest angle between two vectors.
  It exists because sometimes the glm::dot product of two normalized vectors
  can end up outside of [-1, 1] (probably due to floating point arithmetic),
  which makes acos return NaN.
  */
  float GetAngle(const vec3 x, const vec3 y)
  {
    float arg = dot(x, y);
    arg = clamp(arg, -1.0f, 1.0f);
    return acos(arg);
  }

  bool NanCheck(const float x) { return x != x; }
}

BasicRod::BasicRod(const float rod_length,
                   const int num_edges,
                   const float radius,
                   const float density)
  : rod_length_(rod_length)
  , num_edges_(num_edges)
  , radius_(radius)
  , density_(density)
  , vertex_mass_(radius_*radius_*pi<float>()*rod_length_*density_/num_vertices_)

  , rest_curvatures_(num_vertices_)

  , vertices_(num_vertices_)
  , prev_vertices_(num_vertices_)
  , velocities_(num_vertices_)
  , edges_(num_edges_)
  , prev_edges_(num_edges_)
  , binorms_(num_vertices_)
  , external_forces_(num_vertices_)
  , internal_forces_(num_vertices_)

  , reference_frames_(num_edges_)
  , material_frames_(num_edges_)

  , constraints_(num_edges_)
  , constr_diagonal_(num_edges_)
  , constr_off_diagonal_(num_edges_-1)
  , reference_twist_(num_vertices_, 0)
  , material_frame_angles_(num_edges_, 0)
  , angle_velocities_(num_edges_, 0)
  , energy_angle_gradient_(num_edges_)
  , energy_angle_hessian_diag_(num_edges_)
  , energy_angle_hessian_off_diag_(num_edges_-1)
  , rest_stretch_(num_edges_, 0)
{
  Init();
}

void BasicRod::Move()
{
  ComputeEdges();
  ComputeBinorms();
  ComputeMaterialFrames();
  ComputeInternalForces();
  IntegrateForces();
  UpdateReferenceFrames();
  ResetExternalForces();

  if (viscous_)
    SetCurrentShapeToRestShape();
}

void BasicRod::IntegrateForces()
{
  for (int i = 0; i < num_vertices_; ++i)
  {
    velocities_[i] +=
      timestep_*vertex_mass_inv_*
      (external_forces_[i] + internal_forces_[i] - damping_*velocities_[i]) +
      timestep_*gravity_;
    vertices_[i] += timestep_*velocities_[i];
  }
  for (int i = 0; i < num_edges_; ++i)
  {
    angle_velocities_[i] -= timestep_*
      (GetTwistEnergyAngleDerivative(i) + GetBendEnergyAngleDerivative(i))/GetCrossSectionalInertia(i) +
      0.15f*angle_velocities_[i];
    material_frame_angles_[i] += timestep_*angle_velocities_[i];
  }
}

void BasicRod::ComputeEdges()
{
  for (int i = 0; i < num_edges_; ++i)
  {
    edges_[i] = vertices_[i+1] - vertices_[i];
  }
}

void BasicRod::EnforceLengthConstraints()
{
  int count = 0;

  ComputeEdges();
  ComputeConstraints();
  ComputeConstraintMatrix();
  prev_vertices_ = vertices_;

  auto minmax = std::minmax_element(constraints_.begin(), constraints_.end());

  while (*minmax.second-*minmax.first > 0.0001f*rod_length_)
  {
    LAPACKE_sptsv(LAPACK_COL_MAJOR,
                  num_edges_,
                  1,
                  &constr_diagonal_[0],
                  &constr_off_diagonal_[0],
                  &constraints_[0],
                  num_edges_);

    std::vector<float>& lambda = constraints_;

    for (int i = 0; i < num_edges_; ++i)
    {
      vertices_[i] -= -vertex_mass_inv_*edge_length_inv_*lambda[i]*edges_[i];
      vertices_[i+1] -= vertex_mass_inv_*edge_length_inv_*lambda[i]*edges_[i];
    }
    ComputeEdges();
    ComputeConstraints();
    ComputeConstraintMatrix();
  }

  UpdateVelocities();
}

void BasicRod::ComputeConstraints()
{
  for (int i = 0; i < num_edges_; ++i)
  {
    constraints_[i] = 0.5f*edge_length_inv_*dot(edges_[i], edges_[i]) - 0.5f*edge_length_;
  }
}

void BasicRod::ComputeConstraintMatrix()
{
  for (int i = 0; i < num_edges_; ++i)
  {
    constr_diagonal_[i] = 2*vertex_mass_inv_*dot(edges_[i], edges_[i])*edge_length_inv_squared_;
  }

  for (int i = 0; i < num_edges_-1; ++i)
  {
    constr_off_diagonal_[i] = -vertex_mass_inv_*dot(edges_[i], edges_[i+1])*edge_length_inv_squared_;
  }
}

void BasicRod::UpdateVelocities()
{
  for (int i = 0; i < num_vertices_; ++i)
  {
    velocities_[i] += (vertices_[i] - prev_vertices_[i])/timestep_;
  }
}

void BasicRod::Init()
{
  for (int i = 0; i < num_vertices_; ++i)
  {
    vertices_[i] = vec3{0.0f, 0.0f, i*edge_length_};
    velocities_[i] = vec3{0.0f,0.0f,0.0f};
    external_forces_[i] = vec3{0.0f, 0.0f, 0.0f};
    internal_forces_[i] = vec3{0.0f, 0.0f, 0.0f};
    binorms_[i] = vec3{0.0f, 0.0f, 0.0f};
    rest_curvatures_[i] = vec2{0.0f, 0.0f};
    //rest_curvatures_[i] = rotate(vec2(0,0.4f), i*2*3.14f/kNumVertices);
  }
  ComputeEdges();
  prev_edges_ = edges_;

  for (int i = 0; i < num_edges_; ++i)
  {
    //rest_curvatures_[i] = vec2{0.1f, 0.0f};
    reference_frames_[i].first = vec3{1.0f, 0.0f, 0.0f};
    reference_frames_[i].second = vec3{0.0f, 1.0f, 0.0f};
  }
}

void BasicRod::UpdateReferenceFrames()
{

  for (int i = 0; i < num_edges_; ++i)
  {
    vec3 prev_tangent = normalize(prev_edges_[i]);
    vec3 new_tangent = normalize(edges_[i]);
    float angle = GetAngle(prev_tangent, new_tangent);
    if (angle > angle_tolerance_)
    {
      vec3 rotation_axis = normalize(cross(prev_tangent, new_tangent));
      reference_frames_[i].first =
        rotate(reference_frames_[i].first, angle, rotation_axis);

      reference_frames_[i].second =
        cross(new_tangent, reference_frames_[i].first);

      prev_edges_[i] = edges_[i];
    }
  }
  ComputeBinorms();
  ComputeReferenceTwist();
}

void BasicRod::ComputeReferenceTwist()
{
  for (int i = 1; i < num_edges_; ++i)
  {
    vec3 prev_tangent = normalize(edges_[i-1]);
    vec3 next_tangent = normalize(edges_[i]);
    vec3 d1 = reference_frames_[i-1].first;
    float angle = GetAngle(prev_tangent, next_tangent);

    if (angle > angle_tolerance_)
    {
      vec3 rotation_axis = normalize(cross(prev_tangent, next_tangent));
      d1 = rotate(d1, angle, rotation_axis);
    }

    d1 = rotate(d1, reference_twist_[i], next_tangent);
    float sign = glm::sign(dot(d1, reference_frames_[i].second));

    if (sign == 0.0f)
      sign = 1.0f;

    reference_twist_[i] -= sign*GetAngle(d1, reference_frames_[i].first);
  }
}

void BasicRod::ComputeMaterialFrames()
{
  for (int i = 0; i < num_edges_; ++i)
  {
    float sin_theta = sinf(material_frame_angles_[i]);
    float cos_theta = cosf(material_frame_angles_[i]);

    material_frames_[i].first =
      cos_theta*reference_frames_[i].first + sin_theta*reference_frames_[i].second;
    material_frames_[i].second =
      -sin_theta*reference_frames_[i].first + cos_theta*reference_frames_[i].second;
  }
}

vec2 BasicRod::GetCurvatureAngleJacobian(const int i, const int j) const
{
  vec2 jacobian{0.0f};
  if (j == i-1)
  {
    jacobian.x = -0.5f*dot(binorms_[i],material_frames_[i-1].first);
    jacobian.y = -0.5f*dot(binorms_[i],material_frames_[i-1].second);
  }
  else if (j == i)
  {
    jacobian.x = -0.5f*dot(binorms_[i],material_frames_[i].first);
    jacobian.y = -0.5f*dot(binorms_[i],material_frames_[i].second);
  }
  return jacobian;
}

/* Computes the second derivative of curvature i w.r.t. the
 * material frame angle corresponding to edge j.
 * j = PREV for edge i-1 and j = NEXT for edge i.
*/
vec2 BasicRod::GetCurvatureAngleJacobianJacobian(const int i, const int j) const
{
  vec2 jacobian{0.0f};

  if (j == i-1)
  {
    jacobian.x = 0.5f*dot(-binorms_[i],material_frames_[i-1].second);
    jacobian.y = 0.5f*dot(binorms_[i],material_frames_[i-1].first);
  }
  else if (j == i)
  {
    jacobian.x = 0.5f*dot(-binorms_[i],material_frames_[i].second);
    jacobian.y = 0.5f*dot(binorms_[i],material_frames_[i].first);
  }
  return jacobian;
}

void BasicRod::ComputeBinorms()
{
  for (int i = 1; i < num_edges_; ++i)
  {
    float divisor =
      length(edges_[i-1])*length(edges_[i]) + dot(edges_[i-1], edges_[i]);

    if (divisor)
      binorms_[i] = 2.0f*(cross(edges_[i-1], edges_[i])) / divisor;
    else
      binorms_[i] = vec3{0.0f};
  }
}

void BasicRod::ComputeInternalForces()
{
  for (int i = 0; i < num_vertices_; ++i)
  {
    internal_forces_[i] =
      -(GetBendEnergyVertexGradient(i) + GetTwistEnergyVertexGradient(i) + GetStretchEnergyVertexGradient(i));
  }
}

void BasicRod::SetCurrentShapeToRestShape()
{
  for (int i = 0; i < num_edges_; ++i)
  {
    rest_curvatures_[i] = GetCurvature(i);
    rest_stretch_[i] = GetStretch(i);
  }
}

vec3 BasicRod::GetBendEnergyVertexGradient(const int i) const
{
  vec3 gradient{0.0f};
  gradient =
    GetBendStiffness(i-1)*GetCurvatureVertexJacobian(i-1, i)*(GetCurvature(i-1) - GetRestCurvature(i-1)) +
    GetBendStiffness(i)*GetCurvatureVertexJacobian(i, i)*(GetCurvature(i) - GetRestCurvature(i)) +
    GetBendStiffness(i+1)*GetCurvatureVertexJacobian(i+1, i)*(GetCurvature(i+1) - GetRestCurvature(i+1));

  return edge_length_inv_*gradient;
}

/* Compute the TRANSPOSED jacobian of curvature i w.r.t. a vertex. */
mat2x3 BasicRod::GetCurvatureVertexJacobian(const int i, const int j) const
{
  return GetCurvatureEdgeJacobian(i, j-1) - GetCurvatureEdgeJacobian(i, j);
}

/* Compute the TRANSPOSED jacobian of curvature i w.r.t. an edge. */
mat2x3 BasicRod::GetCurvatureEdgeJacobian(const int i, const int j) const
{
  mat2x3 jacobian;

  if (i < 1 || i > num_edges_-1)
    return jacobian;

  vec3 prev_tangent = normalize(edges_[i-1]);
  vec3 next_tangent = normalize(edges_[i]);
  float chi = 1.0f + dot(prev_tangent, next_tangent);

  if (chi == 0.0f)
    return jacobian;

  vec3 tangent_tilde = (prev_tangent + next_tangent)/chi;
  std::pair<vec3, vec3> frame_tilde =
    {(material_frames_[i-1].first + material_frames_[i].first)/chi,
     (material_frames_[i-1].second + material_frames_[i].second)/chi};

  if (j == i-1)
  {
    jacobian[0] =
        (-GetCurvature(i).x*tangent_tilde +
          cross(next_tangent, frame_tilde.second)) / length(edges_[i-1]);
    jacobian[1] =
        (-GetCurvature(i).y*tangent_tilde -
          cross(next_tangent, frame_tilde.first)) / length(edges_[i-1]);
  }
  else if (j == i)
  {
    jacobian[0] =
        (-GetCurvature(i).x*tangent_tilde -
          cross(prev_tangent, frame_tilde.second)) / length(edges_[i]);
    jacobian[1] =
        (-GetCurvature(i).y*tangent_tilde +
          cross(prev_tangent, frame_tilde.first)) / length(edges_[i]);
  }
  else
    jacobian = mat2x3{0.0f};

  return jacobian;
}

vec3 BasicRod::GetTwistEnergyVertexGradient(const int i) const
{
  vec3 gradient{0.0f};

  gradient = GetTwistStiffness(i-1)*GetTwist(i-1)*GetTwistVertexGradient(i-1, i) +
             GetTwistStiffness(i)*GetTwist(i)*GetTwistVertexGradient(i, i) +
             GetTwistStiffness(i+1)*GetTwist(i+1)*GetTwistVertexGradient(i+1, i);

  return edge_length_inv_*gradient;
}

vec3 BasicRod::GetTwistVertexGradient(const int i, const int j) const
{
  return GetTwistEdgeGradient(i, j-1) - GetTwistEdgeGradient(i, j);
}

vec3 BasicRod::GetTwistEdgeGradient(const int i, const int j) const
{
  vec3 gradient{0.0f};

  if (i < 1 || i > num_edges_-1)
    return gradient;

  if (j == i-1)
    gradient = 0.5f*binorms_[i]/length(edges_[i-1]);
  else if (j == i)
    gradient = 0.5f*binorms_[i]/length(edges_[i]);

  return gradient;
}

vec3 BasicRod::GetStretchEnergyVertexGradient(const int i) const
{
  if (i > 0 && i < num_edges_)
  {
    return
      GetStretchStiffness(i-1)*GetStretchVertexGradient(i-1, i)*(GetStretch(i-1) - rest_stretch_[i-1])*edge_length_ +
      GetStretchStiffness(i)*GetStretchVertexGradient(i, i)*(GetStretch(i) - rest_stretch_[i])*edge_length_;
  }
  else if (i == 0)
  {
    return
      GetStretchStiffness(i)*GetStretchVertexGradient(i, i)*(GetStretch(i) - rest_stretch_[i])*edge_length_;
  }
  else if (i == num_edges_)
  {
    return
      GetStretchStiffness(i-1)*GetStretchVertexGradient(i-1, i)*(GetStretch(i-1) - rest_stretch_[i-1])*edge_length_;
  }
  else
    return {0.0f, 0.0f, 0.0f};
}

vec3 BasicRod::GetStretchVertexGradient(const int i, const int j) const
{
  vec3 gradient{0.0f};
  if (j == i)
    gradient = -GetStretchEdgeGradient(i);
  else if (j == i+1)
    gradient = GetStretchEdgeGradient(i);
  return gradient;
}

vec3 BasicRod::GetStretchEdgeGradient(const int i) const
{
  if (i < 0 || i > num_edges_-1)
    return vec3{0.0f};

  return edge_length_inv_*normalize(edges_[i]);
}

float BasicRod::GetStretch(const int i) const
{
  if (i < 0 || i > num_edges_-1)
    return 0.0f;

  return length(edges_[i])*edge_length_inv_ - 1.0f;
}

float BasicRod::GetTwistEnergyAngleDerivative(const int i) const
{
  float derivative =
    GetTwistStiffness(i)*GetTwist(i) - GetTwistStiffness(i+1)*GetTwist(i+1);
  return edge_length_inv_*derivative;
}

float BasicRod::GetBendEnergyAngleDerivative(const int i) const
{
  vec2 curv = GetCurvature(i) - GetRestCurvature(i);
  float derivative =
    GetBendStiffness(i)*dot(curv, GetCurvatureAngleJacobian(i, i));
  curv = GetCurvature(i+1) - GetRestCurvature(i+1);
  derivative +=
    GetBendStiffness(i+1)*dot(curv, GetCurvatureAngleJacobian(i+1, i));

  return edge_length_inv_*derivative;
}

float BasicRod::GetTwist(const int i) const
{
  if (i > 0 && i < num_edges_)
  {
    return
      material_frame_angles_[i] -
      material_frame_angles_[i-1] +
      reference_twist_[i];
  }
  else
    return 0.0f;
}

vec2 BasicRod::GetCurvature(const int i) const
{
  if (i > 0 && i < num_edges_)
  {
    vec2 curvature;
    curvature.x =
      0.5f*(dot(material_frames_[i-1].second + material_frames_[i].second, binorms_[i]));
    curvature.y =
      -0.5f*(dot(material_frames_[i-1].first + material_frames_[i].first, binorms_[i]));
    return curvature;
  }
  else
    return vec2{0.0f};
}

vec2 BasicRod::GetRestCurvature(const int i) const
{
  if (i > 0 && i < num_edges_)
    return rest_curvatures_[i];
  else
    return vec2{0.0f};
}

float BasicRod::GetStretchStiffness(const int i) const
{
  if (viscous_)
    return 3.0f*viscosity_*GetEdgeCrossSectionArea(i)/timestep_;
  else
    return young_modulus_*GetEdgeCrossSectionArea(i);
}

float BasicRod::GetBendStiffness(const int i) const
{
  float radius = GetVertexRadius(i);
  if (viscous_)
    return 3.0f*viscosity_*GetVertexCrossSectionArea(i)*radius*radius/(4.0f*timestep_);
  else
    return 0.25f*young_modulus_*GetVertexCrossSectionArea(i)*radius*radius;
}

float BasicRod::GetTwistStiffness(const int i) const
{
  float radius = GetVertexRadius(i);
  if (viscous_)
    return 0.5f*viscosity_*GetVertexCrossSectionArea(i)*radius*radius/timestep_;
  else
    return 0.5f*shear_modulus_*GetVertexCrossSectionArea(i)*radius*radius;
}

float BasicRod::GetEdgeCrossSectionArea(const int i) const
{
  float length_ratio = edge_length_/length(edges_[i]);
  return length_ratio*pi<float>()*radius_*radius_;
}

float BasicRod::GetVertexCrossSectionArea(const int i) const
{
  float radius = GetVertexRadius(i);
  return pi<float>()*radius*radius;
}

float BasicRod::GetEdgeRadius(const int i) const
{
  if (i < 0 || i > num_edges_-1)
    return 0.0f;
  else
    return radius_*glm::sqrt(edge_length_/length(edges_[i]));
}

float BasicRod::GetVertexRadius(const int i) const
{
  if (i > 0 && i < num_edges_)
    return 0.5f*(GetEdgeRadius(i-1) + GetEdgeRadius(i));
  else if (i == 0)
    return GetEdgeRadius(0);
  else if (i == num_edges_)
    return GetEdgeRadius(num_edges_-1);
  else
    return 0.0f;
}

float BasicRod::GetCrossSectionalInertia(const int i) const
{
  float radius = GetEdgeRadius(i);
  return density_*pi<float>()*radius_*radius_*edge_length_*radius*radius;
}

void BasicRod::Displace(const int i, const glm::vec3& d)
{
  vertices_[i] += d;
  velocities_[i] += d/timestep_;
}

void BasicRod::ResetExternalForces()
{
  for (int i = 0; i < num_vertices_; ++i)
    external_forces_[i] = vec3{0.0f};
}

mat4 BasicRod::GetMaterialFrame(const int i) const
{
  mat4 frame{1.0f};

  frame[0] = vec4{material_frames_[i].first, 0.0f};
  frame[1] = vec4{material_frames_[i].second, 0.0f};
  frame[2] = vec4{normalize(edges_[i]), 0.0f};
  return frame;
}
