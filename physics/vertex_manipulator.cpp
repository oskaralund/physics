#include "vertex_manipulator.h"
#include "vertex_selector.h"
#include "oskgl.h"
#include "distance_point_ray.h"
#include <iostream>

using namespace glm;

VertexManipulator::VertexManipulator(DynamicalSystemInterface* ds,
                                     CameraInterface* cam)
  : dynamical_system_(ds)
  , camera_(cam)
{
  vertex_selector_ = std::make_unique<VertexSelector>(dynamical_system_);
}

void VertexManipulator::Update()
{
  if (selected_vertex_index_ == -1)
    return;

  vec3 selected_vertex = dynamical_system_->GetVertex(selected_vertex_index_);
  dynamical_system_->AddForce(selected_vertex_index_,
                              picking_force_*(target_ - selected_vertex));
  //dynamical_system_->Displace(selected_vertex_index_,
  //                            target_ - selected_vertex);
}

void VertexManipulator::MouseMoveCallback(double x, double y)
{
  if (selected_vertex_index_ == -1)
    return;

  mat2x3 ray = GetWindowRay(x, y);

  vec3 selected_vertex = dynamical_system_->GetVertex(selected_vertex_index_);
  DistancePointRay(selected_vertex, ray[0], ray[1], &target_);
}

void VertexManipulator::MouseDownCallback(double x, double y)
{
  mat2x3 ray = GetWindowRay(x, y);
  selected_vertex_index_ =
    vertex_selector_->GetVertexClosestToRay(ray[0], ray[1]);
  if (selected_vertex_index_ != -1)
  {
    vec3 selected_vertex = dynamical_system_->GetVertex(selected_vertex_index_);
    DistancePointRay(selected_vertex, ray[0], ray[1], &target_);
  }
}

void VertexManipulator::MouseUpCallback(double x, double y)
{
  selected_vertex_index_ = -1;
}

mat2x3 VertexManipulator::GetWindowRay(double x, double y)
{
  mat4 projection_matrix = camera_->GetProjectionMatrix();
  mat4 view_matrix = camera_->GetViewMatrix();
  return oskgl::CastRay(x, y,
                        projection_matrix,
                        view_matrix);
}
