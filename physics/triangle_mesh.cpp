#include "triangle_mesh.h"

using namespace glm;

void TriangleMesh::AddVertex(vec3 p)
{
  vertices_.push_back(p);
  ++num_vertices_;
  vertex_normals_.resize(num_vertices_);
}

void TriangleMesh::AddTriangle(const int i, const int j, const int k)
{
  triangles_.push_back(uvec3{i, j, k});
  ++num_triangles_;
}

void TriangleMesh::ComputeVertexNormals()
{
  for (int n = 0; n < num_triangles_; ++n)
  {
    int i = triangles_[n][0];
    int j = triangles_[n][1];
    int k = triangles_[n][2];
    vec3 surface_normal = ComputeSurfaceNormal(i, j, k);
    vertex_normals_[i] += surface_normal;
    vertex_normals_[j] += surface_normal;
    vertex_normals_[k] += surface_normal;
  }

  for (int n = 0; n < num_vertices_; ++n)
    vertex_normals_[n] = normalize(vertex_normals_[n]);
}

vec3 TriangleMesh::ComputeSurfaceNormal(const int i, const int j, const int k)
{
  vec3 edge1 = vertices_[j] - vertices_[i];
  vec3 edge2 = vertices_[k] - vertices_[i];
  return normalize(cross(edge1, edge2));
}
