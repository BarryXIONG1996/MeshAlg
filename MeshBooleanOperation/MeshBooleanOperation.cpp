#include <Geometry.h>
#include "TopoTriMesh.h"
#include "BooleanOperation.h"
#include <array>
#include <osg/Geometry>
#include <osg/Geode>
#include <osg/Vec3>
#include <osg/Array>
#include <osg/PrimitiveSet>
#include <osg/Group>
#include <osgViewer/Viewer>
#include <osgGA/TrackballManipulator>

// Helper function to add a cube centered at (cx, cy, cz) with side length 'side'
void AddCube(TriMesh& mesh, double cx, double cy, double cz, double side) {
    // Define the vertices of the cube
    std::vector<Vec3d> vertices = {
        {cx - side / 2, cy - side / 2, cz - side / 2},
        {cx + side / 2, cy - side / 2, cz - side / 2},
        {cx + side / 2, cy + side / 2, cz - side / 2},
        {cx - side / 2, cy + side / 2, cz - side / 2},// 底面四点
        {cx - side / 2, cy - side / 2, cz + side / 2},
        {cx + side / 2, cy - side / 2, cz + side / 2},
        {cx + side / 2, cy + side / 2, cz + side / 2},
        {cx - side / 2, cy + side / 2, cz + side / 2} // 顶面四点
    };

    // Store the current size of points to adjust indices accordingly
    size_t baseIndex = 1;

    // Add vertices to the mesh
    for (const auto& vertex : vertices) {
        mesh.points.push_back(vertex);
    }

    // Define the faces of the cube as triangles
    std::vector<std::array<int, 4>> faces = {
        {0, 3, 2, 1}, // bottom face
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

// 将 TriMesh 转为 osg::Geometry
osg::ref_ptr<osg::Geometry> createGeometryFromTriMesh(const TriMesh& mesh)
{
    if (mesh.points.empty() || mesh.indices.empty()) {
        return nullptr;
    }

    // === 1. 创建顶点数组 ===
    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
    vertices->reserve(mesh.points.size());
    for (const auto& p : mesh.points) {
        vertices->push_back(osg::Vec3(
            static_cast<float>(p.x),
            static_cast<float>(p.y),
            static_cast<float>(p.z)
        ));
    }

    // === 2. 解析索引（1-based + 0 分隔）===
    std::vector<GLuint> triIndices; // 使用 GLuint 而非 unsigned int（更明确）
    triIndices.reserve(mesh.indices.size());

    for (size_t i = 0; i < mesh.indices.size(); ) {
        if (mesh.indices[i] == 0) {
            ++i;
            continue;
        }
        if (i + 2 >= mesh.indices.size()) break;

        // 转为 0-based
        GLint i0 = mesh.indices[i] - 1;
        GLint i1 = mesh.indices[i + 1] - 1;
        GLint i2 = mesh.indices[i + 2] - 1;

        // 边界检查
        if (i0 >= 0 && i1 >= 0 && i2 >= 0 &&
            i0 < static_cast<GLint>(mesh.points.size()) &&
            i1 < static_cast<GLint>(mesh.points.size()) &&
            i2 < static_cast<GLint>(mesh.points.size())) {
            triIndices.push_back(static_cast<GLuint>(i0));
            triIndices.push_back(static_cast<GLuint>(i1));
            triIndices.push_back(static_cast<GLuint>(i2));
        }
        i += 4; // 跳过 3 索引 + 1 个 0
    }

    if (triIndices.empty()) {
        return nullptr; // 或返回空几何体
    }

    // === 3. 创建 DrawElementsUInt ===
    osg::ref_ptr<osg::DrawElementsUInt> drawElements = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES);
    drawElements->reserve(triIndices.size());
    for (GLuint idx : triIndices) {
        drawElements->push_back(idx); // 安全：类型匹配 GLuint
    }

    // === 4. 构建 Geometry ===
    osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry;
    geometry->setVertexArray(vertices.get());
    geometry->addPrimitiveSet(drawElements.get());

    // === 5. （可选）生成法向 ===
    osg::ref_ptr<osg::Vec3Array> normals = new osg::Vec3Array;
    normals->resize(vertices->size(), osg::Vec3(0, 0, 0));

    for (size_t i = 0; i < triIndices.size(); i += 3) {
        GLuint a = triIndices[i], b = triIndices[i + 1], c = triIndices[i + 2];
        if (a >= normals->size() || b >= normals->size() || c >= normals->size()) continue;

        const osg::Vec3& va = (*vertices)[a];
        const osg::Vec3& vb = (*vertices)[b];
        const osg::Vec3& vc = (*vertices)[c];
        osg::Vec3 faceNorm = (vb - va) ^ (vc - va);
        (*normals)[a] += faceNorm;
        (*normals)[b] += faceNorm;
        (*normals)[c] += faceNorm;
    }

    for (osg::Vec3& n : *normals) {
        if (n.length2() > 0.0f) n.normalize();
    }
    geometry->setNormalArray(normals.get());
    geometry->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);

    // === 6. （可选）设置颜色 ===
    osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(0.8f, 0.8f, 1.0f, 1.0f));
    geometry->setColorArray(colors.get());
    geometry->setColorBinding(osg::Geometry::BIND_OVERALL);

    return geometry;
}

// 主函数：创建场景
osg::ref_ptr<osg::Node> createSceneFromTriMesh(const TriMesh& mesh) {
    auto geode = new osg::Geode;
    auto geom = createGeometryFromTriMesh(mesh);
    if (geom) {
        geode->addDrawable(geom);
    }

    // 启用光照（需要法向）
    geode->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::ON);
    geode->getOrCreateStateSet()->setMode(GL_LIGHT0, osg::StateAttribute::ON);

    return geode;
}

void ShowTriMesh(const TriMesh& mesh)
{
    if (mesh.points.empty()) {
        return;
    }

    osg::ref_ptr<osg::Node> root = createSceneFromTriMesh(mesh);
    if (!root) {
        return;
    }

    osgViewer::Viewer viewer;
    viewer.setSceneData(root);
    viewer.setCameraManipulator(new osgGA::TrackballManipulator());
    viewer.realize(); // 可选：提前创建窗口

    viewer.run();
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