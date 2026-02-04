#include <Geometry.h>
#include <GeomCalc.h>
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
#include <osg/LineWidth>
#include <osg/Material>
#include <osg/Depth>

#include <fstream>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <exception>
#include <numeric>

#if 0 // 正确性测试
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

osg::ref_ptr<osg::Geode> createGeometryFromTriMesh(const TriMesh& mesh)
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

    // === 2. 解析三角形索引（1-based + 0 分隔）===
    std::vector<GLuint> triIndices;
    for (size_t i = 0; i < mesh.indices.size(); ) {
        if (mesh.indices[i] == 0) {
            ++i;
            continue;
        }
        if (i + 2 >= mesh.indices.size()) break;

        GLint i0 = mesh.indices[i] - 1;
        GLint i1 = mesh.indices[i + 1] - 1;
        GLint i2 = mesh.indices[i + 2] - 1;

        if (i0 >= 0 && i1 >= 0 && i2 >= 0 &&
            i0 < static_cast<GLint>(mesh.points.size()) &&
            i1 < static_cast<GLint>(mesh.points.size()) &&
            i2 < static_cast<GLint>(mesh.points.size())) {
            triIndices.push_back(static_cast<GLuint>(i0));
            triIndices.push_back(static_cast<GLuint>(i1));
            triIndices.push_back(static_cast<GLuint>(i2));
        }
        i += 4; // 跳过 3 个顶点 + 1 个分隔符（0）
    }

    if (triIndices.empty()) {
        return nullptr;
    }

    // === 3. 创建三角形 PrimitiveSet ===
    osg::ref_ptr<osg::DrawElementsUInt> triangles = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES);
    triangles->assign(triIndices.begin(), triIndices.end());

    // === 4. 创建线框 PrimitiveSet（每条边）===
    osg::ref_ptr<osg::DrawElementsUInt> lines = new osg::DrawElementsUInt(osg::PrimitiveSet::LINES);
    lines->reserve(triIndices.size() * 2); // 每个三角形 3 条边 × 2 顶点
    for (size_t i = 0; i < triIndices.size(); i += 3) {
        GLuint a = triIndices[i], b = triIndices[i + 1], c = triIndices[i + 2];
        lines->push_back(a); lines->push_back(b);
        lines->push_back(b); lines->push_back(c);
        lines->push_back(c); lines->push_back(a);
    }

    // === 5. 创建面 Geometry ===
    osg::ref_ptr<osg::Geometry> faceGeom = new osg::Geometry;
    faceGeom->setVertexArray(vertices.get());
    faceGeom->addPrimitiveSet(triangles.get());

    // 法向
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
    for (auto& n : *normals) {
        if (n.length2() > 0.0f) n.normalize();
    }
    faceGeom->setNormalArray(normals.get());
    faceGeom->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);

    // 面颜色（浅蓝）
    osg::ref_ptr<osg::Vec4Array> faceColor = new osg::Vec4Array;
    faceColor->push_back(osg::Vec4(0.8f, 0.8f, 1.0f, 1.0f));
    faceGeom->setColorArray(faceColor.get());
    faceGeom->setColorBinding(osg::Geometry::BIND_OVERALL);

    // === 6. 创建线框 Geometry ===
    osg::ref_ptr<osg::Geometry> wireGeom = new osg::Geometry;
    wireGeom->setVertexArray(vertices.get());
    wireGeom->addPrimitiveSet(lines.get());

    // 线框颜色（红色）
    osg::ref_ptr<osg::Vec4Array> wireColor = new osg::Vec4Array;
    wireColor->push_back(osg::Vec4(1.0f, 0.0f, 0.0f, 1.0f));
    wireGeom->setColorArray(wireColor.get());
    wireGeom->setColorBinding(osg::Geometry::BIND_OVERALL);

    // 线宽 + 深度设置（关键！）
    osg::StateSet* wireSS = wireGeom->getOrCreateStateSet();
    wireSS->setAttribute(new osg::LineWidth(2.0f), osg::StateAttribute::ON);
    wireSS->setMode(GL_CULL_FACE, osg::StateAttribute::OFF); // 关闭背面剔除

    // 关闭深度写入，避免被三角面遮挡
    osg::ref_ptr<osg::Depth> depth = new osg::Depth;
    depth->setWriteMask(false); // 不写深度缓冲
    wireSS->setAttributeAndModes(depth, osg::StateAttribute::ON);

    // 可选：关闭光照，确保颜色准确显示
    wireSS->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

    // === 7. 组装到 Geode ===
    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    geode->addDrawable(faceGeom.get());
    geode->addDrawable(wireGeom.get());

    return geode;
}
// 主函数：创建场景
osg::ref_ptr<osg::Node> createSceneFromTriMesh(const TriMesh& mesh) {
    auto geode = createGeometryFromTriMesh(mesh);

    // 启用光照（需要法向）
    /*geode->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::ON);
    geode->getOrCreateStateSet()->setMode(GL_LIGHT0, osg::StateAttribute::ON);*/

    return geode;
}

void ShowTriMesh(const std::vector<TriMesh>& meshs)
{
    if (meshs.empty()) {
        return;
    }

    osg::ref_ptr<osg::Group> root = new osg::Group();
    for (auto const& mesh : meshs) {
        if (mesh.indices.empty())
            continue;
        osg::ref_ptr<osg::Node> node = createSceneFromTriMesh(mesh);
        root->addChild(node);
    } 
    if (root->getNumChildren() < 1)
        return;

    osgViewer::Viewer viewer;
    viewer.setSceneData(root);
    viewer.setCameraManipulator(new osgGA::TrackballManipulator());
    viewer.realize(); // 可选：提前创建窗口

    viewer.run();
}

// 假设 Vec3dCmp 已在项目中定义（使用 g_epsilon 容差比较）
// 若未定义，需补充：struct Vec3dCmp { bool operator()(const Vec3d& a, const Vec3d& b) const { ... } };

void ShowTriMeshWithVertexColor(
    const std::vector<std::vector<Vec3d>>& tris,
    const std::vector<Vec3d>& boundary,
    const std::vector<std::vector<Vec3d>>& intSegs)
{
    if (tris.empty()) return;

    // === 1. 构建关键点集合（boundary + intSegs，容差去重）===
    std::set<Vec3d, Vec3dCmp> keyPoints;
    auto addPoints = [&keyPoints](const auto& container) {
        for (const auto& p : container) keyPoints.insert(p);
        };

    addPoints(boundary);
    for (const auto& poly : intSegs) addPoints(poly);

    // === 2. 构建OSG几何数据 ===
    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
    osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
    vertices->reserve(tris.size() * 3);
    colors->reserve(tris.size() * 3);

    const osg::Vec4 yellow(1.0f, 1.0f, 0.0f, 1.0f); // boundary/intSegs 中的点
    const osg::Vec4 blue(0.0f, 0.0f, 1.0f, 1.0f);    // 其他点（含CDT生成的交点）

    for (const auto& tri : tris) {
        if (tri.size() != 3) continue;
        for (const auto& v : tri) {
            vertices->push_back(osg::Vec3(static_cast<float>(v.x),
                static_cast<float>(v.y),
                static_cast<float>(v.z)));
            colors->push_back(keyPoints.count(v) ? yellow : blue);
        }
    }

    if (vertices->empty()) return;

    // === 3. 创建带顶点颜色的几何体 ===
    osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
    geom->setVertexArray(vertices);
    geom->setColorArray(colors);
    geom->setColorBinding(osg::Geometry::BIND_PER_VERTEX);

    // 添加三角形图元（连续顶点）
    geom->addPrimitiveSet(new osg::DrawArrays(GL_TRIANGLES, 0, vertices->size()));

    // 禁用光照确保颜色准确显示
    osg::ref_ptr<osg::StateSet> ss = geom->getOrCreateStateSet();
    ss->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
    ss->setMode(GL_BLEND, osg::StateAttribute::ON);
    ss->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);

    // === 4. 构建场景并显示 ===
    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    geode->addDrawable(geom);

    osg::ref_ptr<osg::Group> root = new osg::Group;
    root->addChild(geode);

    osgViewer::Viewer viewer;
    viewer.setSceneData(root);
    viewer.setCameraManipulator(new osgGA::TrackballManipulator());
    viewer.realize();
    viewer.run();
}

int main() {
    TriMesh mesh1, mesh2;
#if 0 // 平面共面情况调试
    mesh1.points = {
        {-1,-1,0},{-1,1,0},{1,1,0},{1,-1,0}
    };
    mesh1.indices = {  1,2,3, 0, 1,3,4, 0 };

    mesh2.points = {
        {0,-1,0},{0,1,0},{2,1,0},{2,-1,0}
    };
    mesh2.indices = { 1,2,3, 0, 1,3,4, 0 };
#endif

#if 0 // 平面相交情况调试
    mesh1.points = {
        {-1,-1,0},{-1,1,0},{1,1,0},{1,-1,0}
    };
    mesh1.indices = { 1,2,3, 0, 1,3,4, 0 };

    mesh2.points = {
        {-1,0,-1},{-1,0,1},{1,0,1},{1,0,-1}
    };
    mesh2.indices = { 1,2,3, 0, 1,3,4, 0 };
#endif

#if 0 // 平面，立方体相交情况调试
    AddCube(mesh1, 0, 0, 0, 2);
    mesh2.points = {
        {-1,-1,0},{-1,1,0},{1,1,0},{1,-1,0}
    };
    mesh2.indices = { 1,2,3, 0, 1,3,4, 0 };
#endif

#if 0 // 立方体相交情况-非共面
    AddCube(mesh1, 0, 0, 0, 2);
    AddCube(mesh2, 1, 1, 1, 2);
#endif 

#if 0 // 立方体相交情况-有共面
    AddCube(mesh1, 0, 0, 0, 2);
    AddCube(mesh2, 1, 1, 0, 2);
#endif

#if 0 // 隧道内空，左侧电缆槽
    mesh1.BuildFromOBJ("TestSamples/leftCCSec3d-1.obj");
    for (int i = 0; i < mesh1.indices.size(); i+=4) {
        std::swap(mesh1.indices.at(i), mesh1.indices.at(i+2));
    }
    mesh2.BuildFromOBJ("TestSamples/insideCrossSecs-1.obj");
#endif

#if 0 // 隧道内空，右侧电缆槽
    mesh1.BuildFromOBJ("TestSamples/rightCCSec3d-1.obj");
    mesh1.FixNormalsToOutside();
    mesh2.BuildFromOBJ("TestSamples/insideCrossSecs-1.obj");
#endif

#if 1 // 两个球体
    mesh1.CreateSphere({0,0,0}, 1, 0);
    mesh2.CreateSphere({1,1,1}, 1, 0);
#endif

    TopoTriMesh* topo1 = new TopoTriMesh;
    topo1->Build(mesh1);
    TopoTriMesh* topo2 = new TopoTriMesh;
    topo2->Build(mesh2);

#ifdef _DRAW
    ShowTriMesh({ mesh1,mesh2 });
#endif

    BooleanOperation booleanOp(topo1, topo2, BooleanOperation::INTERSECTION);
    TopoTriMesh res;
    booleanOp.Execute(res);
     
    TriMesh rM;
    res.ToMesh(rM);
    bool isSolid = GeomCalc::IsClosedSolid(rM);

    return 0;
}
#else 效率测试
// 辅助函数：统计TriMesh中的三角形数量（基于0分隔符）
inline int CountTriangles(const TriMesh& mesh) {
    return static_cast<int>(std::count(mesh.indices.begin(), mesh.indices.end(), 0));
}

// 辅助函数：将微秒转换为带单位的字符串（自动选择ms/μs）
std::string FormatDuration(double microseconds) {
    if (microseconds >= 1000.0)
        return std::to_string(microseconds / 1000.0) + " ms";
    return std::to_string(microseconds) + " μs";
}

void BenchmarkBooleanIntersection(const std::string& csvPath = "boolean_benchmark.csv",
    int maxLevel = 4,
    int iterations = 3)
{
    std::ofstream csv(csvPath);
    if (!csv.is_open()) {
        std::cerr << "❌ 无法创建CSV文件: " << csvPath << std::endl;
        return;
    }

    // 严格按需求精简CSV列：仅保留必要性能与质量指标
    csv << "input_triangles_per_sphere,output_triangles,avg_duration_ms,"
        << "min_duration_ms,max_duration_ms,is_closed_solid,notes\n";

    std::cout << "\n[Boolean Operation Benchmark]\n";
    std::cout << "输出文件: " << csvPath << "\n";
    std::cout << "测试范围: subdivisionLevel 0 → " << maxLevel << " (iterations=" << iterations << ")\n";
    std::cout << "CSV列说明: 每球输入面数, 输出面数, 平均/最小/最大耗时(ms), 封闭性, 备注\n\n";

    for (int level = 0; level <= maxLevel; ++level) {
        // 生成测试网格
        TriMesh mesh1, mesh2;
        mesh1.CreateSphere({ 0, 0, 0 }, 1.0, level);
        mesh2.CreateSphere({ 1, 1, 1 }, 1.0, level);

        int triPerSphere = CountTriangles(mesh1); // 单球三角形数
        int totalInputTris = triPerSphere * 2;    // 仅用于控制台显示

        bool isClosed = false;
        std::string notes = "";
        std::vector<double> durations;
        TriMesh rM;

        // 多轮测试取统计值
        for (int iter = 0; iter < iterations; ++iter) {
            try {
                TopoTriMesh topo1, topo2;
                topo1.Build(mesh1), topo2.Build(mesh2);
                TopoTriMesh res;
                
                BooleanOperation booleanOp(&topo1, &topo2, BooleanOperation::INTERSECTION);

                auto start = std::chrono::high_resolution_clock::now();
                booleanOp.Execute(res);
                auto end = std::chrono::high_resolution_clock::now();

                // 仅首次验证结果质量（避免重复计算）
                if (iter == 0) {
                    res.ToMesh(rM);
                    isClosed = GeomCalc::IsClosedSolid(rM);
                    if (!isClosed) notes = "WARNING: Result not closed solid";
                }

                double duration_us = std::chrono::duration<double, std::micro>(end - start).count();
                durations.push_back(duration_us);
            }
            catch (const std::exception& e) {
                notes = "ERROR: " + std::string(e.what());
                durations.push_back(-1.0);
                break;
            }
            catch (...) {
                notes = "ERROR: Unknown exception";
                durations.push_back(-1.0);
                break;
            }
        }

        // 计算统计量
        int outputTris = (durations[0] >= 0 && !durations.empty()) ? CountTriangles(rM) : 0;
        double avg_us = -1.0, min_us = -1.0, max_us = -1.0;
        if (!durations.empty() && durations[0] >= 0) {
            avg_us = std::accumulate(durations.begin(), durations.end(), 0.0) / iterations;
            min_us = *std::min_element(durations.begin(), durations.end());
            max_us = *std::max_element(durations.begin(), durations.end());
        }

        // ✅ 严格按需求输出CSV：移除 subdivision_level 和 total_input_triangles
        csv << triPerSphere << ","
            << outputTris << ","
            << std::fixed << std::setprecision(6) << (avg_us / 1000.0) << ","
            << std::fixed << std::setprecision(6) << (min_us / 1000.0) << ","
            << std::fixed << std::setprecision(6) << (max_us / 1000.0) << ","
            << (isClosed ? "true" : "false") << ","
            << notes << "\n";

        // 控制台保留完整进度信息（含level和总量，便于人工观察）
        std::cout << "[Level " << std::setw(2) << level << "] "
            << "Input/Sphere: " << std::setw(5) << triPerSphere
            << " | Total Input: " << std::setw(5) << totalInputTris
            << " | Output: " << std::setw(5) << outputTris
            << " | Avg: " << FormatDuration(avg_us)
            << (notes.empty() ? "" : " | " + notes) << "\n";
    }

    csv.close();
    std::cout << "\n✅ 基准测试完成！结果已保存至: " << csvPath << "\n";
    std::cout << "💡 分析建议: 用Excel绘制 'input_triangles_per_sphere vs avg_duration_ms' 曲线\n";
}

// ============ 使用示例 ============

int main() {
    // 标准测试：细分级别0～4，每级3次迭代
    BenchmarkBooleanIntersection("boolean_perf.csv", 6, 6);
    return 0;
}
#endif