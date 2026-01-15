#include "GeomCalc.h"
#include "clipper2/clipper.h"
#include "CDT/include/CDT.h"
#include <cmath>
#include <memory>

Vec3d GeomCalc::CompuateNormal(std::vector<Vec3d> const& pnts)
{
    if (pnts.size() < 3) return Vec3d{ 0, 0, 0 };

    Vec3d normal{ 0, 0, 0 };
    size_t n = pnts.size();
    for (size_t i = 0; i < n; ++i) {
        const auto& curr = pnts[i];
        const auto& next = pnts[(i + 1) % n];
        normal.x += (curr.y - next.y) * (curr.z + next.z);
        normal.y += (curr.z - next.z) * (curr.x + next.x);
        normal.z += (curr.x - next.x) * (curr.y + next.y);
    }

    double len = normal.Length();
    if (len < g_epsilon) {
        return Vec3d{ 0, 0, 0 }; // 退化（共线或重合）
    }

    return normal.Normalization();
}

bool GeomCalc::IsPointOnSegment(const Vec3d& p, const Vec3d& a, const Vec3d& b, double& paramP)
{
    const double eps = g_epsilon;
    const double epsSq = eps * eps;

    Vec3d ab = b - a;
    Vec3d ap = p - a;

    double lenABsq = ab.LengthSq();

    // 情况1：A 和 B 重合（退化线段）
    if (lenABsq < epsSq) {
        // 此时线段就是一个点，判断 p 是否等于 a
        if (ap.LengthSq() < epsSq) {
            paramP = 0.0;
            return true;
        }
        return false;
    }

    // 计算投影参数 t = (ap · ab) / |ab|^2
    double t = ap.Dot(ab) / lenABsq;

    // 将 t 限制在 [0, 1] 范围内（避免因浮点误差略微超出）
    // 但我们仍要检查原始 t 是否在合理范围内
    if (t < -eps || t > 1.0 + eps) {
        return false;
    }

    // 计算投影点：a + t * ab
    Vec3d proj = a + ab * t;
    Vec3d diff = p - proj;

    // 检查 p 到线段的垂直距离是否足够小
    if (diff.LengthSq() > epsSq) {
        return false;
    }

    // 成功！将 paramP 限制在 [0, 1] 内（可选，但更安全）
    paramP = std::max(0.0, std::min(1.0, t));
    return true;
}
double GeomCalc::Point2PlaneDistatnce(Vec3d const& pnt, Vec3d const& o, Vec3d const& n)
{
    return std::fabs(Point2PlaneSignDistance(pnt, o, n));
}

double GeomCalc::Point2PlaneSignDistance(Vec3d const& pnt, Vec3d const& o, Vec3d const& n)
{
    return n.Dot(Vec3d(pnt - o));
}

bool GeomCalc::CalPlanePlaneIntersection(
    Vec3d const& o1, Vec3d const& n1,
    Vec3d const& o2, Vec3d const& n2,
    Vec3d& intO,
    Vec3d& intDir)
{
    const double eps = g_epsilon;
    Vec3d cross = n1.Cross(n2);
    double denom = cross.Dot(cross);

    if (denom < eps) {
        // 平面平行
        return false; // 无论共面与否，按你原逻辑不处理
    }

    intDir = cross / std::sqrt(denom); // 单位方向

    double d1 = n1.Dot(o1);
    double d2 = n2.Dot(o2);
    Vec3d temp = n1 * d2 - n2 * d1;
    intO = cross.Cross(temp) / denom;

    return true;
}

bool GeomCalc::CalLineSegmentIntersection(
    Vec3d const& origin,      // 无限直线起点
    Vec3d const& direction,   // 无限直线方向
    Vec3d const& pt1,         // 线段起点
    Vec3d const& pt2,         // 线段终点
    double& lnP1,             // 输出：直线参数 t
    double& SegP2)            // 输出：线段参数 u ∈ [0,1]
{
    const double eps = g_epsilon;

    // --- Step 1: 检查退化 ---
    Vec3d segDir = pt2 - pt1;
    double segLen2 = segDir.Dot(segDir);
    double lineLen2 = direction.Dot(direction);

    if (segLen2 < eps || lineLen2 < eps) {
        return false; // 线段或直线退化为点
    }

    // --- Step 2: 判断两直线是否平行 ---
    Vec3d cross = direction.Cross(segDir);
    if (cross.Length() < eps) {
        return false; // 平行（含共线），按题意不处理
    }

    // --- Step 3: 构造辅助平面 ---
    // 平面过 'origin'，法向为 cross = direction × segDir
    // 此平面包含原直线（因 direction perp cross）
    Vec3d planeOrigin = origin;
    Vec3d planeNormal = direction.Cross(cross); // 关键：用叉积作为平面法向

    // --- Step 4: 求线段所在直线与该平面的交点参数 u ---
    double u;
    if (!linePlaneIntersect(pt1, segDir, planeOrigin, planeNormal, u)) {
        // 理论上不会发生，因为 segDir · cross = |direction × segDir|^2 != 0
        return false;
    }

    // --- Step 5: 检查交点是否在线段上 ---
    if (u < -eps || u > 1.0 + eps) {
        return false;
    }

    // --- Step 6: 验证交点到原直线的距离（可选但推荐）---
    Vec3d intersectPnt = pt1 + segDir * u;
    Vec3d w = intersectPnt - origin;
    double dist = direction.Cross(w).Length() / std::sqrt(lineLen2);
    if (dist > eps) {
        return false;
    }

    // --- Step 7: 计算直线参数 t ---
    double t = w.Dot(direction) / lineLen2;

    lnP1 = t;
    SegP2 = std::clamp(u, 0.0, 1.0);
    return true;
}

// 计算直线与平面的交点参数
// 直线: P(t) = lineOrigin + t * lineDir
// 平面: (X - planeOrigin) · planeNormal = 0
// 返回: 是否相交；若相交，t 为交点参数
bool GeomCalc::linePlaneIntersect(
    const Vec3d& lineOrigin,
    const Vec3d& lineDir,
    const Vec3d& planeOrigin,
    const Vec3d& planeNormal,
    double& t)
{
    const double eps = g_epsilon;
    double denom = lineDir.Dot(planeNormal);
    if (std::abs(denom) < eps) {
        return false; // 直线与平面平行
    }
    t = (planeOrigin - lineOrigin).Dot(planeNormal) / denom;
    return true;
}

using namespace Clipper2Lib;

bool GeomCalc::TriRegionSplit(
    const std::vector<Vec3d>& tri1,
    const std::vector<Vec3d>& tri2,
    std::vector<std::vector<Vec3d>>& tri1OutTri2,   // tri1 \ tri2
    std::vector<std::vector<Vec3d>>& tri2OutTri1,   // tri2 \ tri1
    std::vector<Vec3d>& tri1CollapseTri2            // tri1 collapse tri2
)
{
    const double eps = g_epsilon;
    const double epsSq = eps * eps;

    // === 1. 输入验证 ===
    if (tri1.size() != 3 || tri2.size() != 3) {
        return false;
    }

    // === 2. 检查退化三角形 ===
    auto isDegenerate = [epsSq](const std::vector<Vec3d>& t) -> bool {
        Vec3d e1 = t[1] - t[0];
        Vec3d e2 = t[2] - t[0];
        Vec3d n = e1.Cross(e2);
        return n.LengthSq() < epsSq;
        };

    if (isDegenerate(tri1) || isDegenerate(tri2)) {
        return false;
    }

    // === 3. 构建共面局部坐标系（以 tri1 为基准）===
    Vec3d origin = tri1[0];
    Vec3d e1 = tri1[1] - tri1[0];
    Vec3d e2 = tri1[2] - tri1[0];
    Vec3d normal = e1.Cross(e2).Normalization();

    // 选择不共线向量构建正交基
    Vec3d arbitrary = (std::abs(normal.x) < 0.9) ? Vec3d{ 1, 0, 0 } : Vec3d{ 0, 1, 0 };
    Vec3d x_axis = normal.Cross(arbitrary).Normalization();
    Vec3d y_axis = normal.Cross(x_axis).Normalization(); // 右手系

    // === 4. 定义 3D to 2D 转换函数 ===
    auto to2D = [&](const Vec3d& p) -> PointD {
        Vec3d rel = p - origin;
        return PointD(rel.Dot(x_axis), rel.Dot(y_axis));
        };

    auto to3D = [&](const PointD& pt) -> Vec3d {
        return origin + x_axis * pt.x + y_axis * pt.y;
        };

    // === 5. 投影三角形到 2D 并确保 CCW ===
    auto projectTri = [&](const std::vector<Vec3d>& tri) -> PathD {
        PathD path = { to2D(tri[0]), to2D(tri[1]), to2D(tri[2]) };
        // 计算有向面积判断方向
        double area = (path[1].x - path[0].x) * (path[2].y - path[0].y) -
            (path[2].x - path[0].x) * (path[1].y - path[0].y);
        if (area < 0) std::reverse(path.begin(), path.end()); // 转为 CCW
        return path;
        };

    PathD p1 = projectTri(tri1);
    PathD p2 = projectTri(tri2);

    PathsD subj{ p1 };
    PathsD clip{ p2 };

    // === 6. 执行布尔运算 ===
    PathsD inter = Intersect(subj, clip, FillRule::EvenOdd);
    PathsD diff1 = Difference(subj, clip, FillRule::EvenOdd); // tri1 \ tri2
    PathsD diff2 = Difference(clip, subj, FillRule::EvenOdd); // tri2 \ tri1

    // === 7. 转换回 3D ===
    auto pathsTo3D = [&](const PathsD& paths) -> std::vector<std::vector<Vec3d>> {
        std::vector<std::vector<Vec3d>> result;
        for (const auto& poly : paths) {
            if (poly.size() < 3) continue; // 忽略退化多边形
            std::vector<Vec3d> poly3D;
            poly3D.reserve(poly.size());
            for (const auto& pt : poly) {
                poly3D.push_back(to3D(pt));
            }
            result.emplace_back(std::move(poly3D));
        }
        return result;
        };

    tri1OutTri2 = pathsTo3D(diff1);
    tri2OutTri1 = pathsTo3D(diff2);

    // === 8. 处理交集（取第一个非空多边形，三角形交集最多一个连通区域）===
    tri1CollapseTri2.clear();
    for (const auto& poly : inter) {
        if (poly.size() >= 3) {
            tri1CollapseTri2.reserve(poly.size());
            for (const auto& pt : poly) {
                tri1CollapseTri2.push_back(to3D(pt));
            }
            break; // 只取第一个（理论上唯一）
        }
    }

    return true;
}

std::vector<std::vector<Vec3d>> GeomCalc::Triangulate(const std::vector<Vec3d>& inputPnts)
{
    if (inputPnts.size() < 3) return {};

    // === 1. 去除首尾重复（闭合多边形常见）===
    std::vector<Vec3d> pnts = inputPnts;
    while (pnts.size() > 2 && (pnts.front() - pnts.back()).LengthSq() <= g_epsilon * g_epsilon)
        pnts.pop_back();
    if (pnts.size() < 3) return {};

    // === 2. 计算法向量（Newell）===
    Vec3d n{ 0, 0, 0 };
    for (size_t i = 0; i < pnts.size(); ++i) {
        const auto& a = pnts[i];
        const auto& b = pnts[(i + 1) % pnts.size()];
        n.x += (a.y - b.y) * (a.z + b.z);
        n.y += (a.z - b.z) * (a.x + b.x);
        n.z += (a.x - b.x) * (a.y + b.y);
    }
    if (n.LengthSq() < g_epsilon * g_epsilon) return {}; // 共线

    // === 3. 选择投影平面 ===
    auto proj = [&](const Vec3d& p) -> CDT::V2d<double> {
        double ax = std::fabs(n.x), ay = std::fabs(n.y), az = std::fabs(n.z);
        if (ax >= ay && ax >= az) return { p.y, p.z }; // YZ
        if (ay >= az)              return { p.x, p.z }; // XZ
        return { p.x, p.y };                           // XY
        };

    // === 4. 构建 2D 点和边界边 ===
    std::vector<CDT::V2d<double>> verts2D;
    verts2D.reserve(pnts.size());
    for (const auto& p : pnts) verts2D.push_back(proj(p));

    // 边：(0,1), (1,2), ..., (n-1,0)
    std::vector<CDT::Edge> edges;
    const size_t N = pnts.size();
    edges.reserve(N);
    for (size_t i = 0; i < N; ++i)
        edges.emplace_back(static_cast<CDT::VertInd>(i), static_cast<CDT::VertInd>((i + 1) % N));

    // === 5. 去重（CDT 工具）===
    CDT::RemoveDuplicatesAndRemapEdges(verts2D, edges);

    if (verts2D.size() < 3) return {};

    // === 6. 执行 CDT ===
    CDT::Triangulation<double> cdt;
    cdt.insertVertices(verts2D);
    cdt.insertEdges(edges);
    cdt.eraseOuterTrianglesAndHoles();

    // === 7. 转回 3D ===
    std::vector<std::vector<Vec3d>> result;
    for (const auto& t : cdt.triangles) {
        // 跳过超级三角形残留（索引越界）
        if (t.vertices[0] >= pnts.size() ||
            t.vertices[1] >= pnts.size() ||
            t.vertices[2] >= pnts.size())
            continue;

        result.push_back({
            pnts[t.vertices[0]],
            pnts[t.vertices[1]],
            pnts[t.vertices[2]]
            });
    }

    return result;
}