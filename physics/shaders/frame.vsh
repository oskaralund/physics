#version 430

layout(std140, binding=0) uniform TransformBlock
{
  mat4 projection;
  mat4 view;
};

layout(location=0) uniform mat4 model;

layout(location=0) in vec4 position;

void main()
{
  gl_Position = projection*view*model*position;
}
