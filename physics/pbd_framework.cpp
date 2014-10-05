#include "pbd_framework.h"

void PBDFramework::CreateVertexBuffers()
{
  const int grid_size = 10;
  glGenBuffers(2, vertices_);
}
