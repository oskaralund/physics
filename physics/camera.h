#ifndef CAMERA_H_
#define CAMERA_H_

#include <glm/glm.hpp>

class Camera
{
public:
  Camera();
  
  void SetPosition(const glm::vec3 &position);
  void SetDirection(const glm::vec3 &direction);
  void SetUp(const glm::vec3 &up);

  void set_fov(const float &fov)             { fov_ = fov; }
  void set_aspect_ratio(const float &aspect) { aspect_ratio_ = aspect; }
  void set_z_near(const float &z_near)       { z_near_ = z_near; }
  void set_z_far(const float &z_far)         { z_far_ = z_far; }

  //Accessors
  glm::mat4 view_matrix() { return view_matrix_; }
  float fov()             { return fov_; }
  float aspect_ratio()    { return aspect_ratio_; }
  float z_near()          { return z_near_; }
  float z_far()           { return z_far_; }

private:
  float fov_          = 0.5f*3.1415f;
  float aspect_ratio_ = 1.0f;
  float z_near_       = 1.0f;
  float z_far_        = 10.0f;
  glm::vec3 position_;
  glm::vec3 direction_;
  glm::vec3 up_;
  glm::mat4 view_matrix_;
  glm::mat4 PV_;
};

#endif
