#version 430

layout(std140, binding=0) uniform TransformBlock
{
  mat4 projection;
  mat4 view;
};

layout(location=0) in vec4 position;
layout(location=1) in vec3 normal;

layout(location=0) uniform mat4 model;
layout(location=1) uniform vec3 dir_to_light;
layout(location=2) uniform vec4 light_intensity;
layout(location=3) uniform mat3 normal_model_to_camera;

smooth out vec4 color;

void main()
{
  gl_Position = projection*view*model*position;
  float cos_ang_incidence = dot(normal, dir_to_light);
  cos_ang_incidence = abs(cos_ang_incidence);
  //color = vec4(1.0,0.0,0.0,1.0);
  color = light_intensity * vec4(1.0, 0.0, 0.0, 1.0) * cos_ang_incidence;
}
