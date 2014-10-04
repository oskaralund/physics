#ifndef VERTEX_MANIPULATOR_INTERFACE_H_
#define VERTEX_MANIPULATOR_INTERFACE_H_

class VertexManipulatorInterface
{
public:
  virtual void Update() = 0; //To be called whenever the dynamical system moves

  /* x, y represent normalized mouse coordinates */
  virtual void MouseMoveCallback(double x, double y) = 0;
  virtual void MouseDownCallback(double x, double y) = 0;
  virtual void MouseUpCallback(double x, double y) = 0;
};

#endif
