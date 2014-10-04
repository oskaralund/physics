#ifndef TRIANGLE_MESH_H_
#define TRIANGLE_MESH_H_

#include <glm/glm.hpp>
#include <vector>
#include <tuple>

class TriangleMesh
{
public:
  void AddVertex(glm::vec3);
  void AddTriangle(const int, const int, const int);
  void ComputeVertexNormals();
  glm::vec3 GetVertex(const int i) const { return vertices_[i]; }
  int GetNumVertices() const { return vertices_.size(); }
  int GetNumTriangles() const { return triangles_.size(); }
  const std::vector<glm::vec3>& GetVertices() { return vertices_; }
  const std::vector<glm::uvec3>& GetTriangles() const { return triangles_; }
  const std::vector<glm::vec3>& GetVertexNormals() const { return vertex_normals_; }
  std::vector<glm::vec3> vertices_;
private:
  int num_vertices_ = 0;
  int num_triangles_ = 0;
  glm::vec3 ComputeSurfaceNormal(const int, const int, const int);
  std::vector<glm::uvec3> triangles_;
  std::vector<glm::vec3> vertex_normals_;
};

#endif
