#ifndef ROD_SANDBOX_H_
#define ROD_SANDBOX_H_

#include "window.h"
#include "basic_rod.h"
#include "rod_vis.h"
#include "fps_camera.h"
#include "grid_vis.h"
#include "vertex_manipulator.h"

#include <memory>

class RodSandbox : public Window
{
public:
  RodSandbox();

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
  BasicRod rod_{0.2f, 20, 0.0025f, 1360.0f};
  RodVis rod_vis_{&rod_};
  float floor_level_ = -0.2f;
  GridVis grid_vis_{1000, 0.1f, glm::vec3{0.0f, -0.2f, 0.0f}, glm::vec3{0.0f, 1.0f, 0.0f}};
  VertexManipulator vertex_manipulator_{&rod_, &camera_};

  GLuint matrix_buffer_;
};

#endif
