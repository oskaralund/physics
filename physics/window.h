#ifndef WINDOW_H_
#define WINDOW_H_

#include <glload/gl_4_3.h>
#include <glload/gl_load.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include <string>
#include <memory>

struct GLFWwindow;

class Window
{
public:
  Window();
  virtual ~Window() {}
  void OpenWindow();
  void SetWindow(Window* i);

protected:
  GLFWwindow* GetGLFWWindowPointer() const { return window_; }
  void AbsoluteToNormalizedCoordinates(double* x, double* y) const;
  void GetNormalizedMouseCoordinates(double* x, double* y) const;
  void CenterCursor();
  void HideCursor(const bool);

  virtual void DrawScene() = 0;
  virtual void UpdateDynamics(double dt) = 0;
  virtual void HandleMouseMove(double x, double y) {} //Assumes normalized mouse coordinates
  virtual void HandleMouseButton(int button, int action, int mods) {}
  virtual void HandleKey(int key, int scancode, int action, int mods) {}
  virtual void Reshape(int width, int height) {}

private:
  int InitializeGLFW();
  inline static void StaticHandleMouseMove(GLFWwindow* w, double x, double y)
  {
    Window* window =
      static_cast<Window*>(glfwGetWindowUserPointer(w));
    window->HandleMouseMove(x, y);
  }
  inline static void StaticHandleMouseButton(GLFWwindow* w, int button, int action, int mods)
  {
    Window* window =
      static_cast<Window*>(glfwGetWindowUserPointer(w));
    window->HandleMouseButton(button, action, mods);
  }
  inline static void StaticHandleKey(GLFWwindow* w, int key, int scancode, int action, int mods)
  {
    Window* window =
      static_cast<Window*>(glfwGetWindowUserPointer(w));
    window->HandleKey(key, scancode, action, mods);
  }
  inline static void StaticReshape(GLFWwindow* w, int width, int height)
  {
    Window* window =
      static_cast<Window*>(glfwGetWindowUserPointer(w));
    window->Reshape(width, height);
  }

  int width_ = 1280;
  int height_ = 720;
  double prev_time_ = 0;
  std::string title_ = "Sandbox";
  GLFWwindow* window_;
};

#endif
