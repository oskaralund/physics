#include "pbd_framework.h"
#include <glm/gtx/norm.hpp>
#include <iostream>
#include <algorithm>

using namespace glm;

namespace
{
  mat3 Skew(vec3 x)
  {
    mat3 skew;
    skew[0] = {0.0f, x[2], -x[1]};
    skew[1] = {-x[2], 0.0f, x[0]};
    skew[2] = {x[1], -x[0], 0.0f};
    return skew;
  }
}

PBDFramework::PBDFramework()
  : velocities_(4, vec3{0.0f})
  , external_forces_(4, vec3{0.0f})
  , inverse_masses_(4, 1.0f)
  , new_vertices_(4)
{
  vertices_ = {{0.0f, 0.0f, 0.0f},
               {0.0f, 0.0f, 1.0f},
               {1.0f, 0.0f, 0.0f},
               {1.0f, 0.0f, 1.0f}};
  num_vertices_ = 4;

  AddLengthConstraint(0, 1, 1.0f);
  AddLengthConstraint(0, 2, 1.0f);
  AddLengthConstraint(1, 3, 1.0f);
  AddLengthConstraint(2, 3, 1.0f);
  AddLengthConstraint(0, 3, glm::sqrt(2.0f));
}

PBDFramework::PBDFramework(const TriangleMesh& mesh)
  : num_vertices_(mesh.GetNumVertices())
  , velocities_(num_vertices_, vec3{0.0f})
  , external_forces_(num_vertices_, vec3{0.0f})
  , inverse_masses_(num_vertices_, num_vertices_/mass_)
  , vertices_(num_vertices_)
  , new_vertices_(num_vertices_)
{
  LoadTriangleMesh(mesh);
}

void PBDFramework::Move()
{
  IntegrateExternalForces();
  DampVelocities();

  for (int i = 0; i < num_vertices_; ++i)
    new_vertices_[i] = vertices_[i] + timestep_*velocities_[i];

  ProjectConstraints();

  UpdateVelocities();
  ResetExternalForces();

  vertices_ = new_vertices_;
}

void PBDFramework::ProjectConstraints()
{
  for (int i = 0; i < num_solver_iterations_; ++i)
  {
    if (num_length_constraints_)
      ProjectLengthConstraints();
    if (num_bend_constraints_)
      ProjectBendConstraints();
  }
}



void PBDFramework::Displace(const int i, const glm::vec3& d)
{
  vertices_[i] += d;
  velocities_[i] += d/timestep_;
}

void PBDFramework::AddForce(const int i, const glm::vec3& f)
{
  external_forces_[i] += f;
}

/* LENGTH CONSTRAINT */
void PBDFramework::ProjectLengthConstraints()
{
  for (int k = 0; k < num_length_constraints_; ++k)
  {
    int i = c_length_dep_vars_[k].first;
    int j = c_length_dep_vars_[k].second;
    const vec3 edge_normal = normalize(new_vertices_[i] - new_vertices_[j]);
    const float length_deviation = length(new_vertices_[i] - new_vertices_[j]) - c_length_lengths_[k];
    const float weight_denominator = inverse_masses_[i] + inverse_masses_[j];
    const float weight1 =
      inverse_masses_[i]/weight_denominator;
    const float weight2 =
      inverse_masses_[j]/weight_denominator;
    float k_prime = 1.0f - pow(1.0f - c_length_stiffness_, 1.0f/num_solver_iterations_);
    new_vertices_[i] -= k_prime*weight1*length_deviation*edge_normal;
    new_vertices_[j] += k_prime*weight2*length_deviation*edge_normal;
  }
}

void PBDFramework::AddLengthConstraint(const int i, const int j, const float l)
{
  ++num_length_constraints_;
  c_length_lengths_.push_back(l);
  c_length_dep_vars_.push_back({i, j});
}
/* END LENGTH CONSTRAINT */

/* BEND CONSTRAINT */
void PBDFramework::ProjectBendConstraints()
{
  for (int k = 0; k < num_bend_constraints_; ++k)
  {
    const auto& dep_vars = c_bend_dep_vars_[k];
    vec3& p1 = new_vertices_[dep_vars[0]];
    vec3& p2 = new_vertices_[dep_vars[1]];
    vec3& p3 = new_vertices_[dep_vars[2]];
    vec3& p4 = new_vertices_[dep_vars[3]];
    const vec3 n1 =
      normalize(cross(p2-p1,p3-p1));
    const vec3 n2 =
      normalize(cross(p2-p1,p4-p1));
    const float d = clamp(dot(n1, n2), -1.0f, 1.0f);

    const vec3 q2 = -(cross(p3-p1, n2) + d*cross(n1, p3-p1))/length(cross(p2-p1, p3-p1))
                    -(cross(p4-p1, n1) + d*cross(n2, p4-p1))/length(cross(p2-p1, p4-p1));
    const vec3 q3 = (cross(p2-p1, n2) + d*cross(n1, p2-p1))/length(cross(p2-p1, p3-p1));
    const vec3 q4 = (cross(p2-p1, n1) + d*cross(n2, p2-p1))/length(cross(p2-p1, p4-p1));
    const vec3 q1 = -q2-q3-q4;

    const float w1 = inverse_masses_[dep_vars[0]];
    const float w2 = inverse_masses_[dep_vars[1]];
    const float w3 = inverse_masses_[dep_vars[2]];
    const float w4 = inverse_masses_[dep_vars[3]];
    const float b = c_bend_angles_[k];
    const float denominator = w1*length2(q1) +
                              w2*length2(q2) +
                              w3*length2(q3) +
                              w4*length2(q4);

    const float epsilon = 1e-4f;
    if (denominator < epsilon)
      continue;

    const float k_prime = 1.0f - pow(1.0f - c_bend_stiffness_, 1.0f/num_solver_iterations_);
    p1 -= k_prime*w1*glm::sqrt(1.0f - d*d)*(acos(d) - b)*q1/denominator;
    p2 -= k_prime*w2*glm::sqrt(1.0f - d*d)*(acos(d) - b)*q2/denominator;
    p3 -= k_prime*w3*glm::sqrt(1.0f - d*d)*(acos(d) - b)*q3/denominator;
    p4 -= k_prime*w4*glm::sqrt(1.0f - d*d)*(acos(d) - b)*q4/denominator;
  }
}

void PBDFramework::AddBendConstraint(std::array<int, 4> dep_vars,
                                     const float b)
{
  ++num_bend_constraints_;
  c_bend_angles_.push_back(b);
  c_bend_dep_vars_.push_back(dep_vars);
}
/* END BEND CONSTRAINT */

void PBDFramework::UpdateVelocities()
{
  for (int i = 0; i < num_vertices_; ++i)
  {
    velocities_[i] = (new_vertices_[i] - vertices_[i])/timestep_;
  }
}

void PBDFramework::ResetExternalForces()
{
  for (int i = 0; i < num_vertices_; ++i)
    external_forces_[i] = vec3{0.0f, 0.0f, 0.0f};
}

void PBDFramework::LoadTriangleMesh(const TriangleMesh& mesh)
{
  for (int i = 0; i < mesh.GetNumVertices(); ++i)
    vertices_[i] = mesh.GetVertex(i);

  /* Add length constraints */
  for (const auto& triangle : mesh.GetTriangles())
  {
    int i = triangle[0];
    int j = triangle[1];
    int k = triangle[2];
    float length;

    length = glm::length(vertices_[i] - vertices_[j]);
    AddLengthConstraint(i, j, length);
    length = glm::length(vertices_[i] - vertices_[k]);
    AddLengthConstraint(i, k, length);
    length = glm::length(vertices_[j] - vertices_[k]);
    AddLengthConstraint(j, k, length);
  }

  /* Add bend constraints */
  const int num_triangles = mesh.GetNumTriangles();
  const auto& triangles = mesh.GetTriangles();

  for (int i = 0; i < num_triangles; ++i)
  {
    for (int j = i+1; j < num_triangles; ++j)
    {
      std::array<int, 4> dependent_variables;
      uvec3 tri_i = triangles[i];
      uvec3 tri_j = triangles[j];
      std::sort(&tri_i[0], &tri_i[0]+3);
      std::sort(&tri_j[0], &tri_j[0]+3);
      auto it = std::set_intersection(&tri_i[0], &tri_i[0]+3,
                                      &tri_j[0], &tri_j[0]+3,
                                      dependent_variables.begin());
      if (it == dependent_variables.begin()+2)
      {
        std::set_symmetric_difference(&tri_i[0], &tri_i[0]+3,
                                      &tri_j[0], &tri_j[0]+3,
                                      dependent_variables.begin()+2);
        AddBendConstraint(dependent_variables, 3.1415f);
      }
    }
  }
}

void PBDFramework::IntegrateExternalForces()
{
  for (int i = 0; i < num_vertices_; ++i)
  {
    velocities_[i] +=
      timestep_*inverse_masses_[i]*external_forces_[i] +
      timestep_*gravity_;
  }
}

void PBDFramework::DampVelocities()
{
  const vec3 center_of_mass = GetCenterOfMass();
  const vec3 center_of_mass_velocity = GetCenterOfMassVelocity();
  const vec3 angular_momentum = GetAngularMomentum();
  for (int i = 0; i < num_vertices_; ++i)
  {
    const vec3 dv = center_of_mass_velocity +
                    cross(angular_momentum, vertices_[i] - center_of_mass) -
                    velocities_[i];
    velocities_[i] += damping_*dv;
  }
}

vec3 PBDFramework::GetAngularMomentum()
{
  const vec3 center_of_mass = GetCenterOfMass();
  vec3 L = {0.0f, 0.0f, 0.0f};
  mat3 I{0.0f};
  for (int i = 0; i < num_vertices_; ++i)
  {
    const vec3 r = vertices_[i] - center_of_mass;
    L += cross(r, velocities_[i]/inverse_masses_[i]);
    const mat3 r_tilde = Skew(r);
    I += r_tilde*transpose(r_tilde)/inverse_masses_[i];
  }
  return inverse(I)*L;
}

vec3 PBDFramework::GetCenterOfMass()
{
  vec3 center_of_mass = {0.0f, 0.0f, 0.0f};
  for (int i = 0; i < num_vertices_; ++i)
    center_of_mass += vertices_[i]/inverse_masses_[i];

  center_of_mass /= mass_;
  return center_of_mass;
}

vec3 PBDFramework::GetCenterOfMassVelocity()
{
  vec3 center_of_mass_velocity = {0.0f, 0.0f, 0.0f};
  for (int i = 0; i < num_vertices_; ++i)
    center_of_mass_velocity += velocities_[i]/inverse_masses_[i];

  center_of_mass_velocity /= mass_;
  return center_of_mass_velocity;
}
