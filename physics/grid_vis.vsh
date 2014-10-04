#version 430 core

layout(std140, binding=0) uniform TransformBlock
{
  mat4 projection;
  mat4 view;
};

layout(location=0) in vec4 position;

layout(location=0) uniform mat4 model;
layout(location=1) uniform int grid_size;
layout(location=2) uniform float square_size;

out vec4 color;

void main()
{
  if ( (gl_InstanceID/grid_size) % 2 == 0)
  {
    if (gl_InstanceID % 2 == 0)
      color = vec4(0.2, 0.2, 0.2, 1.0);
    else
      color = vec4(0.4, 0.4, 0.4, 1.0);
  }
  else
  {
    if (gl_InstanceID % 2 == 1)
      color = vec4(0.2, 0.2, 0.2, 1.0);
    else
      color = vec4(0.4, 0.4, 0.4, 1.0);
  }

  vec4 translated_position = position;
  translated_position.x += gl_InstanceID/grid_size - grid_size/2;
  translated_position.z -= gl_InstanceID%grid_size - grid_size/2;
  gl_Position = projection*view*model*translated_position;
}
