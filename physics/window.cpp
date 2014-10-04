#include "window.h"

Window::Window()
{
  InitializeGLFW();
}

void Window::OpenWindow()
{
  glEnable(GL_DEPTH_TEST);
  while (!glfwWindowShouldClose(window_))
  {
    double time = glfwGetTime();
    double dt = time - prev_time_;
    prev_time_ = time;

    UpdateDynamics(dt);

    glClearColor(0.0f, 0.7f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    DrawScene();

    glfwSwapBuffers(window_);
    glfwPollEvents();
  }
  glfwTerminate();
}

int Window::InitializeGLFW()
{
  if (!glfwInit())
    return -1;

  window_ = glfwCreateWindow(width_, height_, title_.c_str(), nullptr, nullptr);

  if (!window_)
  {
    glfwTerminate();
    return -1;
  }

  glfwMakeContextCurrent(window_);
  glfwSetWindowUserPointer(window_, this);

  ogl_LoadFunctions();

  glfwSetMouseButtonCallback(window_, StaticHandleMouseButton);
  glfwSetCursorPosCallback(window_, StaticHandleMouseMove);
  glfwSetKeyCallback(window_, StaticHandleKey);
  glfwSetWindowSizeCallback(window_, StaticReshape);

  return 0;
}

void Window::AbsoluteToNormalizedCoordinates(double* x, double* y) const
{
  int w, h;
  glfwGetWindowSize(window_, &w, &h);
  *x = 2.0*(*x)/w - 1.0;
  *y = 1.0 - 2.0*(*y)/h;
}

void Window::GetNormalizedMouseCoordinates(double* x, double* y) const
{
  int w, h;
  glfwGetWindowSize(window_, &w, &h);
  glfwGetCursorPos(window_, x, y);
  *x = 2.0*(*x)/w - 1.0;
  *y = 1.0 - 2.0*(*y)/h;
}

void Window::CenterCursor()
{
  int w, h;
  glfwGetWindowSize(window_, &w, &h);
  glfwSetCursorPos(window_, w/2, h/2);
}

void Window::HideCursor(const bool hide)
{
  if (hide)
    glfwSetInputMode(window_, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
  else
    glfwSetInputMode(window_, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
}
