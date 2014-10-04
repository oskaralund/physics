#ifndef ROD_VIS_H_
#define ROD_VIS_H_

#include <glm/glm.hpp>
#include "cylinder_vis.h"
#include "frame_vis.h"

class RodInterface;

class RodVis
{
public:
  RodVis(RodInterface* rod) : rod_(rod) {}
  void Draw();
  void DrawMaterialFrames(const bool);

private:
  RodInterface* rod_;
  CylinderVis edge_{100};
  FrameVis frame_{true, true, false};
  bool draw_material_frames_ = false;
};

#endif
