#include <string>
#include <glload/gl_4_3.h>
#include <glm/glm.hpp>

namespace oskgl
{
  char*  ReadFile(const std::string &filename);
  GLuint CompileVertexShader(const std::string& filename);
  GLuint CompileFragmentShader(const std::string& filename);
  bool   GetCompileStatus(const GLint shader);
  GLuint CompileShaderProgram(const std::string& vshader, const std::string& fshader);

  /* Returns a ray in the form of a glm::mat2x3 object, where column 0 is the ray origin
   * and column 1 is the normalized direction.
  */
  glm::mat2x3 CastRay(double mouse_x,
                      double mouse_y,
                      const glm::mat4& projection,
                      const glm::mat4& view);
}
