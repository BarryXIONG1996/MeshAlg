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

class RTreeAccelerator : public SpatialAccelerator {
public:
    virtual void Build(std::vector<Face*> const& faces) override {
        // 构建R树
        for (Face* face : faces) {
            // 将每个面片的包围盒插入到R树中
            m_rtree.Insert((double*)&face->bbox.lowerBnd, (double*)&face->bbox.upperBnd, face);
        }
    }

    virtual std::vector<Face*> Query(const BndBox3d& box) override {
        // 查询R树
        std::vector<Face*> result;
        m_rtree.Search((double*)&box.lowerBnd, (double*)&box.upperBnd,
            [&result](Face* face)->bool {
                result.push_back(face);
                return true; // 继续搜索
            });
        return result;
    }

    virtual void Clear() override {
        m_rtree.RemoveAll();
    }

    virtual void Remove(Face* f) override {
        m_rtree.Remove((double*)&f->bbox.lowerBnd, (double*)&f->bbox.upperBnd, f);
    }

private:
    RTree<Face*, double, 3> m_rtree;
};

MeshIntersector::MeshIntersector(TopoTriMesh& m1, TopoTriMesh& m2)
    : m_mesh1(m1), m_mesh2(m2)
{
    m_accelerator = std::make_shared<RTreeAccelerator>();
}

// 判断点 p 是否在 Edge e 上（不包括端点）
static bool IsPointOnEdgeInterior(const Vec3d& p, Edge* e)
{
    if (!e || !e->v1 || !e->v2) return false;

    const Vec3d& a = e->v1->pnt;
    const Vec3d& b = e->v2->pnt;

    const double eps = g_epsilon;

    // Step 1: 检查点是否与端点重合（排除端点）
    if ((p - a).Length() < eps || (p - b).Length() < eps) {
        return false; // 在端点上，不算“内部”
    }

    // Step 2: 检查三点是否共线：(b - a) × (p - a) ≈ 0
    Vec3d ab = b - a;
    Vec3d ap = p - a;
    Vec3d cross = ab.Cross(ap);
    if (cross.Length() > eps * ab.Length()) {
        return false; // 不共线
    }

    // Step 3: 检查点是否在线段内部（开区间）
    double len2_ab = ab.Dot(ab);
    if (len2_ab < eps) {
        return false; // 边退化为点（不应发生）
    }

    // 计算投影参数 t ∈ [0,1]：p = a + t * ab
    double t = ap.Dot(ab) / len2_ab;

    // 开区间：0 < t < 1
    return (t > eps && t < 1.0 - eps);
}

// 按线段顺序输出边上点（去重 + 排序）
static std::vector<Vec3d> SortPointsOnEdge(
    const Vec3d& A,
    const Vec3d& B,
    const std::vector<Vec3d>& pointsOnEdge)
{
    std::vector<Vec3d> result;
    if (pointsOnEdge.empty()) return result;

    Vec3d AB = B - A;
    double len2_AB = AB.Dot(AB);

    // 处理退化线段（A == B）
    if (len2_AB < g_epsilon * g_epsilon) {
        return { A }; // 所有点都重合于 A
    }

    // 计算每个点的 t 参数
    struct PointWithT {
        Vec3d p;
        double t;
        bool operator<(const PointWithT& other) const {
            return t < other.t;
        }
    };
    std::vector<PointWithT> ptsWithT;
    ptsWithT.reserve(pointsOnEdge.size());
    for (const auto& P : pointsOnEdge) {
        Vec3d AP = P - A;
        double t = AP.Dot(AB) / len2_AB;
        // 可选：钳制 t 到 [0,1]（若点略微超出）
        // t = std::max(0.0, std::min(1.0, t));
        ptsWithT.push_back({ P, t });
    }

    // 按 t 排序
    std::sort(ptsWithT.begin(), ptsWithT.end());

    // 去重（基于 t 的 epsilon 比较）
    result.push_back(ptsWithT.front().p);
    for (size_t i = 1; i < ptsWithT.size(); ++i) {
        if (ptsWithT[i].t - ptsWithT[i - 1].t > g_epsilon) {
            result.push_back(ptsWithT[i].p);
        }
    }

    return result;
}

bool MeshIntersector::Execute(TopoTriMesh& coPlanes)
{
    auto IntersectProcess = [&](std::function<void(Face*, Face*)> intersectFunc, bool coPlanar) {
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
        m_accelerator->Clear();
        m_accelerator->Build(intersectingFaces);

        // 2、相交测试
        for (int idx = 0; idx < m_mesh1.fs.size(); ++idx)
        {
            Face* face = m_mesh1.fs.at(idx);
            // 过滤掉与包围盒不相交的三角面片
            if (face->bbox.IsOut(intBox))
                continue;
            // 查询空间搜索树，获取可能相交的面片
            std::vector<Face*> candidateFaces = m_accelerator->Query(face->bbox);
            for (Face* candidateFace : candidateFaces) {
                // 3、面片之间的相交计算
                FaceFaceInt(face, candidateFace, intersectFunc, coPlanar);
            }
        }
    };

    // 首先, 计算共面部分并移除
    IntersectProcess(std::bind(&MeshIntersector::CoPlanarFaceInt, this, std::placeholders::_1, std::placeholders::_2), true);
    coPlanes = m_coPlanes;

    // 然后, 非共面部分进行求交
    IntersectProcess(std::bind(&MeshIntersector::NonCoPlanarFaceInt, this, std::placeholders::_1, std::placeholders::_2), false);
    
    // 最后，基于求交结果编辑拓扑
    std::set<Face*> rmFaces;
    std::set<Edge*> rmEdges;
    for (auto const& [f, segs] : m_face2Segs) {
        auto const& es = f->getEdges();
        for (auto seg : segs) {
            int segS = m_intSegs.at(seg).first;
            int segE = m_intSegs.at(seg).second;
            for (auto const& e : es) {
                if (IsPointOnEdgeInterior(m_intersectPnts.at(segS), e)) 
                    m_edge2Ints[e].insert(segS);
                if (IsPointOnEdgeInterior(m_intersectPnts.at(segE), e)) 
                    m_edge2Ints[e].insert(segE);
            }
        }
    }

    for (auto const& [f, _] : m_face2Segs) {
        rmFaces.insert(f);
        std::vector<Edge*> es = f->getEdges();
        for (auto* e : es) {
            if (m_edge2Ints.count(e)) {
                if (e->lF == f) rmFaces.insert(e->rF);
                else rmFaces.insert(e->lF);
            }
        }
    }
    for (auto const& [e, _] : m_edge2Ints) rmEdges.insert(e);

    std::vector<std::vector<Vec3d>> reTris1, reTris2;
    for (auto const& f : rmFaces) {
        auto const& es = f->getEdges();
        std::vector<Vec3d> boundary;
        for (auto const& e : es) {
            Vec3d start,end;
            if (e->lF == f)
                start = e->v1->pnt, end = e->v2->pnt;
            else
                start = e->v2->pnt, end = e->v1->pnt;
            std::vector<Vec3d> interPnts;
            if (m_edge2Ints.count(e))
                std::for_each(m_edge2Ints.at(e).begin(), m_edge2Ints.at(e).end(), [&](int idx) { interPnts.push_back(m_intersectPnts.at(idx)); });
            boundary.push_back(start);
            boundary.insert(boundary.end(), interPnts.begin(), interPnts.end());
            boundary.push_back(end); // 首尾点重复添加，之后会去重
        }
        std::vector<std::vector<Vec3d>> intSegs;
        for (auto const& seg : m_face2Segs.at(f))
        {
            Vec3d s = m_intersectPnts.at(m_intSegs.at(seg).first);
            Vec3d e = m_intersectPnts.at(m_intSegs.at(seg).second);
            intSegs.push_back({ s,e });
        }
        std::vector<std::vector<Vec3d>> tris = GeomCalc::TriangulateWithConstraints(boundary, intSegs);
        if (f->topo == &m_mesh1)
            reTris1.insert(reTris1.end(), tris.begin(), tris.end());
        if (f->topo == &m_mesh2)
            reTris2.insert(reTris2.end(), tris.begin(), tris.end());
    }

    for (auto& f : rmFaces) {
        if (f->topo == &m_mesh1) m_mesh1.RemoveFace(f);
        else m_mesh2.RemoveFace(f);
        if (f) delete f;
    }
    for (auto& e : rmEdges) {
        if (e->lF) {
            if (e->lF->topo == &m_mesh1){
                m_mesh1.RemoveEdge(e);
            } else {
                m_mesh2.RemoveEdge(e);
            }
        }
        else {
            if (e->rF->topo == &m_mesh1) {
                m_mesh1.RemoveEdge(e);
            }
            else {
                m_mesh2.RemoveEdge(e);
            }
        }
        if (e) delete e;
    }

    for (auto& tri : reTris1) m_mesh1.AddFace2TopoTriMesh(tri);
    for (auto& tri : reTris2) m_mesh2.AddFace2TopoTriMesh(tri);

    return true;
}

std::vector<Vec3d> MeshIntersector::EdgeEdgeInt(Edge* e1, Edge* e2)
{
    return std::vector<Vec3d>();
}

void MeshIntersector::FaceFaceInt(Face* f1, Face* f2, std::function<void(Face*, Face*)> intersectFunc, bool coPlanar)
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
        if (coPlanar)
            intersectFunc(f1, f2);
    }
    else
    {
        // 面不平行
        if (!coPlanar)
            intersectFunc(f1, f2);
    }

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

    for (auto& e : e2Pts1) {
        Face* cF = (e.first->lF == f1) ? e.first->rF : e.first->lF;
        if (cF) {
            std::vector<Vec3d> newPoly = BuildPolygonWithEdgePoints(cF, e.first, e.second);
            poly1s.push_back(std::move(newPoly));
        }
        m_mesh1.RemoveFace(cF);
    }
    m_mesh1.RemoveFace(f1);
    for (auto& e : e2Pts1) m_mesh1.RemoveEdge(e.first);


    for (auto& e : e2Pts2) {
        Face* cF = (e.first->lF == f2) ? e.first->rF : e.first->lF;
        if (cF) {
            std::vector<Vec3d> newPoly = BuildPolygonWithEdgePoints(cF, e.first, e.second);
            poly2s.push_back(std::move(newPoly));
        }
        m_mesh2.RemoveFace(cF);
        m_accelerator->Remove(cF);
    }
    m_mesh2.RemoveFace(f2);
    m_accelerator->Remove(f2);
    for (auto& e : e2Pts2) m_mesh2.RemoveEdge(e.first);

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
        {
            Face* newF = m_mesh2.AddFace2TopoTriMesh(tri);
            m_accelerator->Build({newF});
        }
    }

    std::vector<std::vector<Vec3d>> tris = GeomCalc::Triangulate(tri12Int);
    for (auto const& tri : tris)
        m_coPlanes.AddFace2TopoTriMesh(tri);
}

enum TopoType { EdgeType, VertexType };

struct IntersectionInfo {
    double param;
    void* topo;
    TopoType type;
};

// 判断顶点和面的相对位置关系
static void  VertexPlaneRelPos(Vertex* v, Vec3d const& o, Vec3d const& n) {
    if (!v) return;
    double signDist = GeomCalc::Point2PlaneSignDistance(v->pnt, o, n);
    if (signDist > g_epsilon) v->posTag = 2; // 0=default, 1=in, 2=out, 3=on
    else if (signDist < -g_epsilon) v->posTag = 1;
    else v->posTag = 3;
};

// 根据两个交点信息构建 efs 和权重
static std::pair<int, std::set<std::pair<Edge*, Face*>>> BuildEfsAndWeight(
    const IntersectionInfo& info1, 
    const IntersectionInfo& info2, 
    Vec3d const& o1, Vec3d const& n1,
    Vec3d const& o2, Vec3d const& n2
)
{
    std::set<std::pair<Edge*, Face*>> efs;
    int weight = 0;

    if (!info1.topo || !info2.topo) return {};

    if (info1.type == TopoType::EdgeType && info2.type == TopoType::EdgeType) {
        weight = 2;
        Edge* e1 = static_cast<Edge*>(info1.topo);
        Edge* e2 = static_cast<Edge*>(info2.topo);
        if (e1->lF) efs.insert({ e2, e1->lF });
        if (e1->rF) efs.insert({ e2, e1->rF });
        if (e2->lF) efs.insert({ e1, e2->lF });
        if (e2->rF) efs.insert({ e1, e2->rF });
        VertexPlaneRelPos(e1->v1, o2, n2);
        VertexPlaneRelPos(e1->v2, o2, n2);
        VertexPlaneRelPos(e2->v1, o1, n1);
        VertexPlaneRelPos(e2->v2, o1, n1);
    }
    else if (info1.type == TopoType::EdgeType && info2.type == TopoType::VertexType) {
        weight = 3;
        Edge* e1 = static_cast<Edge*>(info1.topo);
        Vertex* v2 = static_cast<Vertex*>(info2.topo);
        for (Edge* ce : v2->GetAdjacentEdges()) {
            if (!ce) continue;
            if (e1->lF) efs.insert({ ce, e1->lF });
            if (e1->rF) efs.insert({ ce, e1->rF });
        }
        for (Face* cf : v2->GetAdjacentFaces()) {
            if (cf) efs.insert({ e1, cf });
        }
        v2->posTag = 3;
        VertexPlaneRelPos(e1->v1, o2, n2);
        VertexPlaneRelPos(e1->v2, o2, n2);
    }
    else if (info1.type == TopoType::VertexType && info2.type == TopoType::EdgeType) {
        weight = 3;
        Vertex* v1 = static_cast<Vertex*>(info1.topo);
        Edge* e2 = static_cast<Edge*>(info2.topo);
        for (Edge* ce : v1->GetAdjacentEdges()) {
            if (!ce) continue;
            if (e2->lF) efs.insert({ ce, e2->lF });
            if (e2->rF) efs.insert({ ce, e2->rF });
        }
        for (Face* cf : v1->GetAdjacentFaces()) {
            if (cf) efs.insert({ e2, cf });
        }
        v1->posTag = 3;
        VertexPlaneRelPos(e2->v1, o1, n1);
        VertexPlaneRelPos(e2->v2, o1, n1);
    }
    else { // Vertex-Vertex
        weight = 4;
        Vertex* v1 = static_cast<Vertex*>(info1.topo);
        Vertex* v2 = static_cast<Vertex*>(info2.topo);
        auto cfs1 = v1->GetAdjacentFaces();
        auto cfs2 = v2->GetAdjacentFaces();
        for (Edge* ce : v1->GetAdjacentEdges()) {
            if (!ce) continue;
            for (Face* cf : cfs2) {
                if (cf) efs.insert({ ce, cf });
            }
        }
        for (Edge* ce : v2->GetAdjacentEdges()) {
            if (!ce) continue;
            for (Face* cf : cfs1) {
                if (cf) efs.insert({ ce, cf });
            }
        }
        v1->posTag = v2->posTag = 3;
    }

    return { weight, efs };
}

// 处理“面”与拓扑元素相交的情况（即一个端点来自另一面内部）
static std::pair<int, std::set<std::pair<Edge*, Face*>>> BuildEfsForFaceIntersect(
    const IntersectionInfo& info, 
    Face* faceOnOtherSide, 
    Vec3d const& oOtherSide, 
    Vec3d const& nOtherSide
)
{
    std::set<std::pair<Edge*, Face*>> efs;
    int weight = (info.type == TopoType::EdgeType) ? 1 : 2;

    if (info.type == TopoType::EdgeType) {
        Edge* e = static_cast<Edge*>(info.topo);
        efs.insert({ e, faceOnOtherSide });
        VertexPlaneRelPos(e->v1, oOtherSide, nOtherSide);
        VertexPlaneRelPos(e->v2, oOtherSide, nOtherSide);
    }
    else {
        Vertex* v = static_cast<Vertex*>(info.topo);
        for (Edge* ce : v->GetAdjacentEdges()) {
            if (ce) {
                efs.insert({ ce, faceOnOtherSide });
                /*VertexPlaneRelPos(ce->v1, oOtherSide, nOtherSide);
                VertexPlaneRelPos(ce->v2, oOtherSide, nOtherSide);*/
            }
        }
        v->posTag = 3;
    }

    return { weight, efs };
}

void MeshIntersector::NonCoPlanarFaceInt(Face* f1, Face* f2)
{
    // --- 平面求交 ---
    auto ps1 = f1->getPnts(), ps2 = f2->getPnts();
    if (ps1.size() < 3 || ps2.size() < 3) return;

    Vec3d o1 = ps1[0], norm1 = GeomCalc::CompuateNormal(ps1);
    Vec3d o2 = ps2[0], norm2 = GeomCalc::CompuateNormal(ps2);
    Vec3d intO, intDir;
    if (!GeomCalc::CalPlanePlaneIntersection(o1, norm1, o2, norm2, intO, intDir))
        return;

    // --- 收集交点参数 ---
    struct DbCompare {
        bool operator()(double a, double b) const { return a < b - g_epsilon; }
    };
    using ParamMap = std::map<double, std::tuple<void*, TopoType>, DbCompare>;
    ParamMap param1, param2;

    auto collectParams = [&](const std::vector<Edge*>& edges, ParamMap& out) {
        for (Edge* e : edges) {
            if (!e) continue;
            auto pnts = e->getPnts();
            if (pnts.size() != 2) continue;
            double LineT, segT;
            if (!GeomCalc::CalLineSegmentIntersection(intO, intDir, pnts[0], pnts[1], LineT, segT))
                continue;
            if (segT < g_epsilon) {
                out[LineT] = { e->v1, TopoType::VertexType };
            }
            else if (segT > 1 - g_epsilon) {
                out[LineT] = { e->v2, TopoType::VertexType };
            }
            else {
                out[LineT] = { e, TopoType::EdgeType };
            }
        }
    };

    collectParams(f1->getEdges(), param1);
    collectParams(f2->getEdges(), param2);

    // --- 有效性检查 ---
    if (param1.size() < 2 || param2.size() < 2) return;
    double t1_min = param1.begin()->first, t1_max = param1.rbegin()->first;
    double t2_min = param2.begin()->first, t2_max = param2.rbegin()->first;
    if (t1_max <= t2_min + g_epsilon || t2_max <= t1_min + g_epsilon) return;

    // --- 提取端点信息 ---
    auto makeInfo = [](const ParamMap::value_type& kv) {
        return IntersectionInfo{ kv.first, std::get<0>(kv.second), std::get<1>(kv.second) };
    };

    IntersectionInfo p1_start = makeInfo(*param1.begin());
    IntersectionInfo p1_end = makeInfo(*param1.rbegin());
    IntersectionInfo p2_start = makeInfo(*param2.begin());
    IntersectionInfo p2_end = makeInfo(*param2.rbegin());

    std::vector<Vec3d> intPnts;
    std::vector<int> weights;
    std::vector<std::set<std::pair<Edge*, Face*>>> allEfs;

    // --- 处理起点 ---
    if (std::abs(p1_start.param - p2_start.param) <= g_epsilon) {
        auto [w, efs] = BuildEfsAndWeight(p1_start, p2_start, o1, norm1, o2, norm2);
        weights.push_back(w);
        allEfs.push_back(efs);
        intPnts.push_back(intO + intDir * p1_start.param);
    }
    else if (p1_start.param < p2_start.param) {
        auto [w, efs] = BuildEfsForFaceIntersect(p2_start, f1, o1, norm1);
        weights.push_back(w);
        allEfs.push_back(efs);
        intPnts.push_back(intO + intDir * p2_start.param);
    }
    else {
        auto [w, efs] = BuildEfsForFaceIntersect(p1_start, f2, o2, norm2);
        weights.push_back(w);
        allEfs.push_back(efs);
        intPnts.push_back(intO + intDir * p1_start.param);
    }

    // --- 处理终点 ---
    if (std::abs(p1_end.param - p2_end.param) < g_epsilon) {
        auto [w, efs] = BuildEfsAndWeight(p1_end, p2_end, o1, norm1, o2, norm2);
        weights.push_back(w);
        allEfs.push_back(efs);
        intPnts.push_back(intO + intDir * p1_end.param);
    }
    else if (p2_end.param < p1_end.param) {
        auto [w, efs] = BuildEfsForFaceIntersect(p2_end, f1, o1, norm1);
        weights.push_back(w);
        allEfs.push_back(efs);
        intPnts.push_back(intO + intDir * p2_end.param);
    }
    else {
        auto [w, efs] = BuildEfsForFaceIntersect(p1_end, f2, o2, norm2);
        weights.push_back(w);
        allEfs.push_back(efs);
        intPnts.push_back(intO + intDir * p1_end.param);
    }

    // --- 存储结果 ---
    if (intPnts.size() != weights.size() || intPnts.size() != allEfs.size() || intPnts.size() != 2)
        return;
    std::vector<int> seg;
    for (int intIdx = 0; intIdx < intPnts.size(); ++intIdx) // 解决奇异性问题
    {
        auto const& efs = allEfs.at(intIdx);
        bool exist = false;
        for (auto const& ef : efs) 
        {
            if (!m_ef2Int.count(ef))
                continue;
            exist = true;
            assert(m_ef2Int.at(ef).size() < 2, "线面交点多于2个，错误！");
            double w = weights.at(intIdx);
            int otherIntIdx = *m_ef2Int.at(ef).begin();
            double otherIntW = m_weights.at(otherIntIdx);
            if (w >= 2.0 && otherIntW >= 2.0) // 线面共面，此时可以有两个交点
            {
                m_intersectPnts.push_back(intPnts.at(intIdx));
                m_weights.push_back(w);
                int segEndIdx = m_intersectPnts.size() - 1;
                seg.push_back(segEndIdx);
            }
            else // 取高权重点
            {
                if (w > otherIntW)
                {
                    m_intersectPnts.at(otherIntIdx) = intPnts.at(intIdx);
                    m_weights.at(otherIntIdx) = w;
                }
                seg.push_back(otherIntIdx);
            }
        }

        if (exist)
        {
            for (auto const& ef : efs)
            {
                m_ef2Int[ef].insert(seg.back());
            }
        }
        else
        {
            m_intersectPnts.push_back(intPnts.at(intIdx));
            m_weights.push_back(weights.at(intIdx));
            int segEndIdx = m_intersectPnts.size() - 1;
            seg.push_back(segEndIdx);
            for (auto const& ef : efs)
            {
                m_ef2Int.insert({ ef, {(int)segEndIdx} });
            }
        }
    }
    m_intSegs.push_back({ seg.front(),seg.back() });
    int intSegIndex = m_intSegs.size() - 1;
    m_face2Segs[f1].insert(intSegIndex);
    m_face2Segs[f2].insert(intSegIndex);
}
