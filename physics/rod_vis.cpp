#include "rod_vis.h"
#include "rod_interface.h"

void RodVis::Draw()
{
  for (int i = 0; i < rod_->GetNumEdges(); ++i)
  {
    edge_.SetStartPoint(rod_->GetVertex(i));
    edge_.SetEndPoint(rod_->GetVertex(i+1));
    edge_.SetRadius(rod_->GetEdgeRadius(i));
    edge_.SetLength(rod_->GetEdgeLength(i));

    edge_.SetOrientation(rod_->GetMaterialFrame(i));
    edge_.Draw();

    if (draw_material_frames_)
    {
      frame_.SetOrigo(0.5f*(rod_->GetVertex(i)+rod_->GetVertex(i+1)));
      frame_.SetOrientation(rod_->GetMaterialFrame(i));
      frame_.SetScale(6.0f*rod_->GetEdgeRadius(i));
      frame_.Draw();
    }
  }

}

void RodVis::DrawMaterialFrames(const bool b)
{
  draw_material_frames_ = b;
}
