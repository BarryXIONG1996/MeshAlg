#include "MeshIntersector.h"
#include <GeomCalc.h>
#include <memory>
#include <cmath>
#include <set>

// 定义 M_PI 和 M_E，防止未声明错误
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#include "RTree.h"

extern const double g_epsilon;

class SpatialAccelerator {
public:
    virtual void Build(std::vector<Face*>& faces) = 0;
    virtual std::vector<Face*> Query(const BndBox3d& box) = 0;
};

class RTreeAccelerator : public SpatialAccelerator {
    virtual void Build(std::vector<Face*>& faces) override 
    {
        // 构建R树
        for (Face* face : faces) {
            // 将每个面片的包围盒插入到R树中
            m_rtree.Insert((double*)&face->bbox.lowerBnd, (double*)&face->bbox.upperBnd, face);
        }
    }

    virtual std::vector<Face*> Query(const BndBox3d& box) override 
    {
        // 查询R树
        std::vector<Face*> result;
        m_rtree.Search((double*)&box.lowerBnd, (double*)&box.upperBnd,
            [&result](Face* face)->bool {
                result.push_back(face);
                return true; // 继续搜索
            });
        return result;
    }

    RTree<Face*, double, 3> m_rtree;
};

MeshIntersector::MeshIntersector(TopoTriMesh& m1, TopoTriMesh& m2)
    : m_mesh1(m1), m_mesh2(m2)
{

}

bool MeshIntersector::Execute(TopoTriMesh& coPlanes)
{
    // 1、构建空间搜索树
    // 计算工具集和目标集的包围盒
    BndBox3d boundingBoxO = m_mesh1.GetBndBox();
    BndBox3d boundingBoxT = m_mesh2.GetBndBox();
    // 计算包围盒的交集
    BndBox3d intBox = boundingBoxO.Intersect(boundingBoxT);
    // 遍历工具集中的所有三角面片
    std::vector<Face*>const& facesT = m_mesh2.fs;
    std::vector<Face*> intersectingFaces;
    for (Face* face : facesT) {
        // 检查三角面片的包围盒是否与intBox相交
        if (!face->bbox.IsOut(intBox)) {
            intersectingFaces.push_back(face);
        }
    }
    // 将相交的面片添加到空间搜索树
    std::shared_ptr<SpatialAccelerator> accelerator = std::make_shared<RTreeAccelerator>();
    accelerator->Build(intersectingFaces);

    // 2、相交测试
    std::vector<Face*>const& facesO = m_mesh1.fs;
    for (Face* face : facesO)
    {
        // 过滤掉与包围盒不相交的三角面片
        if (face->bbox.IsOut(intBox))
            continue;
        // 查询空间搜索树，获取可能相交的面片
        std::vector<Face*> candidateFaces = accelerator->Query(face->bbox);
        for (Face* candidateFace : candidateFaces) {
            // 3、面片之间的相交计算
            FaceFaceInt(face, candidateFace);
        }
    }
    
    coPlanes = m_coPlanes;

    return true;
}

std::vector<Vec3d> MeshIntersector::EdgeEdgeInt(Edge* e1, Edge* e2)
{
    return std::vector<Vec3d>();
}

void MeshIntersector::FaceFaceInt(Face* f1, Face* f2)
{
    if (f1->bbox.IsOut(f2->bbox))
        return;

    std::vector<Vec3d> ps1 = f1->getPnts(), ps2 = f2->getPnts();
    if (ps1.size() < 3 || ps2.size() < 3)
        return;

    Vec3d n1 = GeomCalc::CompuateNormal(ps1), n2 = GeomCalc::CompuateNormal(ps2);

    if (n1.Length() < g_epsilon || n2.Length() < g_epsilon)
        return;

    if (n1.Parallel(n2)) // 面平行
    {
        if (GeomCalc::Point2PlaneDistatnce(ps2.at(0), ps1.at(0), n1) > g_epsilon
            || GeomCalc::Point2PlaneDistatnce(ps1.at(0), ps2.at(0), n2) > g_epsilon)
            return; // 不共面

        // 共面
        CoPlanarFaceInt(f1, f2);

        return;
    }

    // 面不平行
    NonCoPlanarFaceInt(f1, f2);

    return;
}

static std::vector<Vec3d> BuildPolygonWithEdgePoints(
    Face* face,
    Edge* edge,
    const std::vector<std::pair<Vec3d, double>>& ptsOnEdge)
{
    if (!face || !edge) return {};

    auto origPts = face->getPnts();     // [v0, v1, v2]
    auto edges = face->getEdges();    // [e0, e1, e2]

    // 找到 edge 在 face 中的索引
    int edgeIdx = -1;
    for (int i = 0; i < 3; ++i) {
        if (edges[i] == edge) {
            edgeIdx = i;
            break;
        }
    }
    if (edgeIdx < 0) return origPts;

    // 面内该边的起点和终点
    const Vec3d& faceFrom = origPts[edgeIdx];
    const Vec3d& faceTo = origPts[(edgeIdx + 1) % 3];

    // 判断面内方向是否与 edge->v1 -> edge->v2 一致
    bool sameDir = (faceFrom.Equal(edge->v1->pnt) && faceTo.Equal(edge->v2->pnt));

    // 按面内顺序排序交点
    auto sorted = ptsOnEdge;
    if (sameDir) {
        std::sort(sorted.begin(), sorted.end(),
            [](const auto& a, const auto& b) { return a.second < b.second; });
    }
    else {
        std::sort(sorted.begin(), sorted.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });
    }

    // 构建新多边形
    std::vector<Vec3d> poly = origPts;
    size_t insertPos = (edgeIdx + 1) % poly.size();
    for (const auto& p : sorted) {
        poly.insert(poly.begin() + insertPos++, p.first);
    }

    return poly;
}

void MeshIntersector::CoPlanarFaceInt(Face* f1, Face* f2)
{
    std::vector<Vec3d> tri1 = f1->getPnts();
    std::vector<Vec3d> tri2 = f2->getPnts();

    std::vector<std::vector<Vec3d>> tri1Out2, tri2Out1;
    std::vector<Vec3d> tri12Int;
    GeomCalc::TriRegionSplit(tri1, tri2, tri1Out2, tri2Out1, tri12Int);

    if (tri12Int.size() < 3) return;

    auto edges1 = f1->getEdges();
    auto edges2 = f2->getEdges();

    std::map<Edge*, std::vector<std::pair<Vec3d, double>>> e2Pts1, e2Pts2;
    for (const Vec3d& pt : tri12Int) {
        for (Edge* e : edges1) {
            double t = 0.0;
            if (GeomCalc::IsPointOnSegment(pt, e->v1->pnt, e->v2->pnt, t) &&
                t > g_epsilon && t < 1 - g_epsilon) {
                e2Pts1[e].push_back({ pt, t });
                break;
            }
        }
        for (Edge* e : edges2) {
            double t = 0.0;
            if (GeomCalc::IsPointOnSegment(pt, e->v1->pnt, e->v2->pnt, t) &&
                t > g_epsilon && t < 1 - g_epsilon) {
                e2Pts2[e].push_back({ pt, t });
                break;
            }
        }
    }

    std::vector<std::vector<Vec3d>> poly1s = std::move(tri1Out2);
    std::vector<std::vector<Vec3d>> poly2s = std::move(tri2Out1);

    m_mesh1.RemoveFace(f1);
    for (auto& e : e2Pts1) {
        Face* cF = (e.first->lF == f1) ? e.first->rF : e.first->lF;
        if (cF) {
            std::vector<Vec3d> newPoly = BuildPolygonWithEdgePoints(cF, e.first, e.second);
            poly1s.push_back(std::move(newPoly));
        }
        m_mesh1.RemoveFace(cF);
        m_mesh1.RemoveEdge(e.first);
    }

    m_mesh2.RemoveFace(f2);
    for (auto& e : e2Pts2) {
        Face* cF = (e.first->lF == f2) ? e.first->rF : e.first->lF;
        if (cF) {
            std::vector<Vec3d> newPoly = BuildPolygonWithEdgePoints(cF, e.first, e.second);
            poly2s.push_back(std::move(newPoly));
        }
        m_mesh2.RemoveFace(cF); 
        m_mesh2.RemoveEdge(e.first);
    }

    for (auto const& poly : poly1s)
    {
        std::vector<std::vector<Vec3d>> tris = GeomCalc::Triangulate(poly);
        for (auto const& tri : tris)
            m_mesh1.AddFace2TopoTriMesh(tri);
    }

    for (auto const& poly : poly2s)
    {
        std::vector<std::vector<Vec3d>> tris = GeomCalc::Triangulate(poly);
        for (auto const& tri : tris)
            m_mesh2.AddFace2TopoTriMesh(tri);
    }

    std::vector<std::vector<Vec3d>> tris = GeomCalc::Triangulate(tri12Int);
    for (auto const& tri : tris)
        m_coPlanes.AddFace2TopoTriMesh(tri);
}

void MeshIntersector::NonCoPlanarFaceInt(Face* f1, Face* f2)
{
}
