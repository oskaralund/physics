#version 430 core

layout(std140, binding=0) uniform TransformBlock
{
  mat4 projection;
  mat4 view;
};

layout(location=0) in vec4 position;
layout(location=1) in vec4 color;

layout(location=0) uniform mat4  model;
layout(location=3) uniform float twist;

out vec4 vs_color;

void main()
{
  mat4 twist_rot = mat4(1.0);
  twist_rot[0].xy = vec2(cos(twist), -sin(twist));
  twist_rot[1].xy = vec2(sin(twist),  cos(twist));
  vec4 pos = position;
  if (position.z == 1)
  {
    pos = twist_rot*pos;
  }
  vs_color    = color;
  gl_Position = projection*view*model*pos;
}
