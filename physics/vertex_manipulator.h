#ifndef VERTEX_MANIPULATOR_H_
#define VERTEX_MANIPULATOR_H_

#include "vertex_manipulator_interface.h"
#include "vertex_selector_interface.h"
#include "dynamical_system_interface.h"
#include "camera_interface.h"
#include <glm/glm.hpp>
#include <memory>

class VertexManipulator : public VertexManipulatorInterface
{
public:
  VertexManipulator(DynamicalSystemInterface* ds,
                    CameraInterface* cam);

  /* Call Update whenever the dynamical system moves */
  void Update() override;

  /* x, y represent normalized mouse coordinates */
  void MouseMoveCallback(double x, double y) override;
  void MouseDownCallback(double x, double y) override;
  void MouseUpCallback(double x, double y) override;

  void SetPickingForce(const float f) { picking_force_ = f; }

private:
  glm::mat2x3 GetWindowRay(double x, double y);

  int selected_vertex_index_ = -1;
  float picking_force_ = 1.0f;
  DynamicalSystemInterface* dynamical_system_;
  CameraInterface* camera_;
  std::unique_ptr<VertexSelectorInterface> vertex_selector_;
  glm::vec3 target_;
};

#endif
