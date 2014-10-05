#ifndef PBD_FRAMEWORK_H_
#define PBD_FRAMEWORK_H_

#include <glload/gl_4_3.h>
#include <glload/gl_load.h>
#include <memory>
#include <array>

class PBDFramework
{
public:

private:
  void CreateVertexBuffers();
  GLuint vertices_[2];
};

#endif
