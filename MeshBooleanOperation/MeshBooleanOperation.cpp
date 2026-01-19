#include <Geometry.h>
#include "TopoTriMesh.h"
#include "BooleanOperation.h"
#include <array>

// Helper function to add a cube centered at (cx, cy, cz) with side length 'side'
void AddCube(TriMesh& mesh, double cx, double cy, double cz, double side) {
    // Define the vertices of the cube
    std::vector<Vec3d> vertices = {
        {cx - side / 2, cy - side / 2, cz - side / 2},
        {cx + side / 2, cy - side / 2, cz - side / 2},
        {cx + side / 2, cy + side / 2, cz - side / 2},
        {cx - side / 2, cy + side / 2, cz - side / 2},
        {cx - side / 2, cy - side / 2, cz + side / 2},
        {cx + side / 2, cy - side / 2, cz + side / 2},
        {cx + side / 2, cy + side / 2, cz + side / 2},
        {cx - side / 2, cy + side / 2, cz + side / 2}
    };

    // Store the current size of points to adjust indices accordingly
    size_t baseIndex = 1;

    // Add vertices to the mesh
    for (const auto& vertex : vertices) {
        mesh.points.push_back(vertex);
    }

    // Define the faces of the cube as triangles
    std::vector<std::array<int, 4>> faces = {
        {0, 1, 2, 3}, // bottom face
        {4, 5, 6, 7}, // top face
        {0, 1, 5, 4}, // front face
        {1, 2, 6, 5}, // right face
        {2, 3, 7, 6}, // back face
        {3, 0, 4, 7}  // left face
    };

    // Add faces to the mesh as triangles
    for (const auto& face : faces) {
        mesh.indices.push_back(baseIndex + face[0]);
        mesh.indices.push_back(baseIndex + face[1]);
        mesh.indices.push_back(baseIndex + face[2]);
        mesh.indices.push_back(0);
        mesh.indices.push_back(baseIndex + face[0]);
        mesh.indices.push_back(baseIndex + face[2]);
        mesh.indices.push_back(baseIndex + face[3]);
        mesh.indices.push_back(0);
    }
}

int main() {
    TriMesh mesh1, mesh2;

    // Add first cube centered at origin with side length 2
    AddCube(mesh1, 0, 0, 0, 2);

    // Add second cube translated along X axis by 1 unit with side length 2
    AddCube(mesh2, 1, 0, 0, 2);

    TopoTriMesh* topo1 = new TopoTriMesh;
    topo1->Build(mesh1);
    TopoTriMesh* topo2 = new TopoTriMesh;
    topo2->Build(mesh2);

    BooleanOperation booleanOp(topo1, topo2, BooleanOperation::DIFFERENCE);
    TopoTriMesh res;
    booleanOp.Execute(res);

    return 0;
}