#define _CRT_SECURE_NO_WARNINGS
#include "oskgl.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <glm/gtx/vector_angle.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace glm;

namespace oskgl {

  GLuint CompileShaderProgram(const std::string& vshader, const std::string& fshader)
  {
    GLuint vshader_id = CompileVertexShader(vshader);
    GLuint fshader_id = CompileFragmentShader(fshader);
    GLuint program_id = glCreateProgram();
    glAttachShader(program_id, vshader_id);
    glAttachShader(program_id, fshader_id);
    glLinkProgram(program_id);
    return program_id;
  }

//TODO: fix this to not use strncpy so we can remove _CRT_SECURE_NO_WARNINGS define.
  char*  ReadFile(const std::string &filename)
  {
    std::ifstream in(filename);
    if (in.is_open())
    {
      std::string contents((std::istreambuf_iterator<char>(in)),
                           (std::istreambuf_iterator<char>()));
      char *sz_contents = new char[contents.size()+1];
      std::strncpy(sz_contents, contents.c_str(), contents.size()+1);
      return sz_contents;
    }
    else
      return nullptr;
  }

  GLuint CompileVertexShader(const std::string& filename)
  {
    char *source = ReadFile(filename);
    if (source == nullptr)
      return -1;

    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, (const GLchar**)&source, NULL);
    glCompileShader(vertex_shader);
    bool compiled_correctly = GetCompileStatus(vertex_shader);

    if (compiled_correctly)
      return vertex_shader;
    else
      return -1;
  }

  GLuint CompileFragmentShader(const std::string &filename)
  {
    char *source = ReadFile(filename);
    if (source == nullptr)
      return -1;

    GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, (const GLchar**)&source, NULL);
    glCompileShader(fragment_shader);
    bool compiled_correctly = GetCompileStatus(fragment_shader);
    if (compiled_correctly)
      return fragment_shader;
    else
      return -1;
  }

  bool GetCompileStatus(const GLint shader)
  {
    GLint compiled;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
    if (compiled)
    {
      return true;
    }
    else
    {
      GLint log_length;
      glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &log_length);
      char* msg_buffer = new char[log_length];
      glGetShaderInfoLog(shader, log_length, NULL, msg_buffer);
      std::cout << msg_buffer << std::endl;
      delete msg_buffer;
      return false;
    }
  }

  mat2x3 CastRay(double x,
                 double y,
                 const mat4& projection,
                 const mat4& view)
  {
    mat2x3 ray;

    vec4 ray_clip{x, y, -1.0f, 1.0f};
    vec4 ray_eye = inverse(projection)*ray_clip;
    ray_eye.z = -1.0f;
    ray_eye.w = 0.0f;
    mat4 camera = inverse(view);
    auto ray_wor = vec3{camera*ray_eye};
    ray[0] = vec3{camera[3]};
    ray[1] = normalize(ray_wor);
    return ray;
  }
}
