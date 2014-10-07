#ifndef FLUID_SANDBOX_H_
#define FLUID_SANDBOX_H_

#include "window.h"
#include "fps_camera.h"
#include "grid_vis.h"

#include <memory>

class FluidSandbox : public Window
{
public:
  FluidSandbox();

protected:
  void DrawScene() override;
  void UpdateDynamics(double elapsed_time) override;
  void HandleKey(int key, int scancode, int action, int mods) override;
  void HandleMouseMove(double x, double y) override;
  void HandleMouseButton(int button, int action, int mods) override;
  void Reshape(int width, int height) override;

private:
  void CreateMatrixBuffer();
  void UpdateViewMatrix();
  void UpdateProjectionMatrix();
  void MoveCamera(const float dt);
  bool move_forward_ = false;
  bool move_back_    = false;
  bool move_left_    = false;
  bool move_right_   = false;
  bool move_up_      = false;
  bool move_down_    = false;
  float camera_sensitivity_ = 0.5f;
  float time_accumulator_ = 0.0f;
  double prev_x_ = 0.0;
  double prev_y_ = 0.0;
  FPSCamera camera_;
  GridVis grid_vis_{1000, 1.0f, glm::vec3{0.0f, -1.0f, 0.0f}, glm::vec3{0.0f, 1.0f, 0.0f}};

  GLuint matrix_buffer_;
};

#endif
