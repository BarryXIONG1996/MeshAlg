#include "GeomCalc.h"
#include "clipper2/clipper.h"
#include "CDT/include/CDT.h"
#include <cmath>
#include <memory>

using namespace Clipper2Lib;

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

// 通过自动选择 XY/XZ/YZ 中面积最大的投影平面
bool GeomCalc::IsLeft(const Vec3d& a, const Vec3d& b, const Vec3d& p, Vec3d const& refN)
{
    Vec3d ab = b - a;
    Vec3d ap = p - a;

    // 计算叉积 (用于找最大投影面)
    Vec3d cross = ab.Cross(ap);

    return cross.Dot(refN) > 0;
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

std::vector<Vec3d> GeomCalc::OrderPointsOnSegment(
    const Vec3d& start,
    const Vec3d& end,
    const std::vector<Vec3d>& pts)
{
    std::vector<Vec3d> result;
    result.reserve(pts.size() + 2);

    // 添加 start 和 end
    result.push_back(start);
    result.push_back(end);
    result.insert(result.end(), pts.begin(), pts.end());

    // 如果 start == end，直接返回去重点（避免除零）
    Vec3d dir = end - start;
    double lenSq = dir.LengthSq();
    if (lenSq < g_epsilon) {
        // 所有点重合，返回唯一代表点
        return { start };
    }

    // 按照从 start 到 end 的方向参数 t 排序
    std::sort(result.begin(), result.end(),
        [&start, &dir, lenSq](const Vec3d& a, const Vec3d& b) {
            // 计算 t = ((p - start) · dir) / |dir|^2
            double t_a = (a - start).Dot(dir) / lenSq;
            double t_b = (b - start).Dot(dir) / lenSq;
            return t_a < t_b;
        });

    // 去重
     std::vector<Vec3d> uniqueResult;
     uniqueResult.push_back(result[0]);
     for (size_t i = 1; i < result.size(); ++i) {
         if ((result[i] - uniqueResult.back()).Length() > g_epsilon) {
             uniqueResult.push_back(result[i]);
         }
     }
     return uniqueResult;
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

// 辅助：计算多边形有向面积（2D）
static double SignedArea(const PathD& p) {
    if (p.size() < 3) return 0.0;
    double area = 0.0;
    size_t n = p.size();
    for (size_t i = 0; i < n; ++i) {
        size_t j = (i + 1) % n;
        area += p[i].x * p[j].y - p[j].x * p[i].y;
    }
    return area * 0.5;
};

// 主函数：多边形求交（共面假设）
bool GeomCalc::PolyIntersect(
    const std::vector<Vec3d>& poly1,
    const std::vector<Vec3d>& poly2,
    std::vector<Vec3d>& outIntersection)
{
    const double eps = g_epsilon;
    const double epsSq = eps * eps;
    outIntersection.clear();

    if (poly1.size() < 3 || poly2.size() < 3) {
        return true; // 无效多边形，无交集
    }

    // === 1. 用 poly1 计算参考法向和平面 ===
    Vec3d origin = poly1[0];
    // 遍历 poly1 的连续三角形扇区，累加法向
    Vec3d refNormal = CompuateNormal(poly1);
    if (refNormal.LengthSq() < epsSq) {
        // poly1 退化（共线或点）
        return true;
    }

    // === 2. 构建局部坐标系 ===
    Vec3d arbitrary = (std::abs(refNormal.x) < 0.9) ? Vec3d{ 1, 0, 0 } : Vec3d{ 0, 1, 0 };
    Vec3d x_axis = refNormal.Cross(arbitrary);
    if (x_axis.LengthSq() < epsSq) {
        arbitrary = Vec3d{ 0, 0, 1 };
        x_axis = refNormal.Cross(arbitrary);
    }
    x_axis = x_axis.Normalization();
    Vec3d y_axis = refNormal.Cross(x_axis); // 自动单位化，右手系

    // === 3. 投影函数 ===
    auto to2D = [&](const Vec3d& p) -> PointD {
        Vec3d rel = p - origin;
        return PointD(rel.Dot(x_axis), rel.Dot(y_axis));
        };

    auto to3D = [&](const PointD& pt) -> Vec3d {
        return origin + x_axis * pt.x + y_axis * pt.y;
        };

    // === 4. 投影并统一为 CCW（相对于 refNormal）===
    auto projectAsCCW = [&](const std::vector<Vec3d>& poly) -> PathD {
        PathD path;
        path.reserve(poly.size());
        for (const auto& p : poly) {
            path.push_back(to2D(p));
        }
        if (path.size() < 3) return path;

        double area = SignedArea(path);
        if (area < 0) {
            std::reverse(path.begin(), path.end());
        }
        return path;
        };

    PathD p1 = projectAsCCW(poly1);
    PathD p2 = projectAsCCW(poly2);

    // === 5. 执行交集运算 ===
    PathsD inter = Intersect({ p1 }, { p2 }, FillRule::EvenOdd);

    // === 6. 转回 3D，取最大面积的多边形（或第一个有效）===
    PathD bestPoly;
    double maxArea = 0.0;

    for (const auto& poly : inter) {
        if (poly.size() >= 3) {
            double area = std::abs(SignedArea(poly));
            if (area > maxArea) {
                maxArea = area;
                bestPoly = poly;
            }
        }
    }

    if (bestPoly.empty()) {
        return true; // 无有效交集
    }

    // 转为 3D
    std::vector<Vec3d> candidate;
    candidate.reserve(bestPoly.size());
    for (const auto& pt : bestPoly) {
        candidate.push_back(to3D(pt));
    }

    // === 7. 确保法向与 poly1 一致 ===
    // 用 candidate 前三个不共线点计算法向
    Vec3d candNormal{ 0, 0, 0 };
    bool foundNormal = false;
    for (size_t i = 0; i + 2 < candidate.size(); ++i) {
        Vec3d e1 = candidate[i + 1] - candidate[i];
        Vec3d e2 = candidate[i + 2] - candidate[i];
        Vec3d n = e1.Cross(e2);
        if (n.LengthSq() > epsSq) {
            candNormal = n;
            foundNormal = true;
            break;
        }
    }

    if (foundNormal && candNormal.Dot(refNormal) < 0) {
        std::reverse(candidate.begin(), candidate.end());
    }

    outIntersection = std::move(candidate);
    return true;
}

    std::vector<std::vector<Vec3d>> GeomCalc::Triangulate(const std::vector<Vec3d>& inputPnts)
    {
        if (inputPnts.size() < 3) return {};

        // === 1. 去除首尾重复 ===
        std::vector<Vec3d> pnts = inputPnts;
        while (pnts.size() > 2 && (pnts.front() - pnts.back()).LengthSq() <= g_epsilon * g_epsilon)
            pnts.pop_back();
        if (pnts.size() < 3) return {};

        // === 2. 计算法向（Newell）===
        Vec3d n = CompuateNormal(pnts);
        if (n.LengthSq() < g_epsilon * g_epsilon) return {}; // 共线

        // === 3. 投影函数 ===
        auto proj = [&](const Vec3d& p) -> CDT::V2d<double> {
            double ax = std::fabs(n.x), ay = std::fabs(n.y), az = std::fabs(n.z);
            if (ax >= ay && ax >= az) return { p.y, p.z }; // YZ
            if (ay >= az)              return { p.x, p.z }; // XZ
            return { p.x, p.y };                           // XY
            };

        // === 4. 构建 2D 点 ===
        std::vector<CDT::V2d<double>> verts2D;
        verts2D.reserve(pnts.size());
        for (const auto& p : pnts) verts2D.push_back(proj(p));

        // === 5. 计算 2D 有符号面积，判断绕向 ===
        size_t N = verts2D.size();
        // 计算标准有符号面积
        double signedArea = 0.0;
        for (size_t i = 0; i < N; ++i) {
            size_t j = (i + 1) % N;
            signedArea += verts2D[i].x * verts2D[j].y - verts2D[j].x * verts2D[i].y;
        }
        bool needReverse = (signedArea < 0); // <0 表示 CW（在标准数学坐标系中）

        // 如果是 CW，反转点序使其变为 CCW（CDT 要求）
        std::vector<Vec3d> pntsForTri = pnts;
        std::vector<CDT::V2d<double>> verts2DForTri = verts2D;
        if (needReverse) {
            std::reverse(pntsForTri.begin(), pntsForTri.end());
            std::reverse(verts2DForTri.begin(), verts2DForTri.end());
        }

        // === 6. 构建边界边（现在是 CCW）===
        std::vector<CDT::Edge> edges;
        edges.reserve(N);
        for (size_t i = 0; i < N; ++i) {
            edges.emplace_back(
                static_cast<CDT::VertInd>(i),
                static_cast<CDT::VertInd>((i + 1) % N)
            );
        }

        // === 7. 去重 ===
        CDT::RemoveDuplicatesAndRemapEdges(verts2DForTri, edges);
        if (verts2DForTri.size() < 3) return {};

        // === 8. CDT 三角化 ===
        CDT::Triangulation<double> cdt;
        cdt.insertVertices(verts2DForTri);
        cdt.insertEdges(edges);
        cdt.eraseOuterTrianglesAndHoles();

        // === 9. 转回 3D，并恢复原始绕向（如果曾反转）===
        std::vector<std::vector<Vec3d>> result;
        for (const auto& t : cdt.triangles) {
            if (t.vertices[0] >= pntsForTri.size() ||
                t.vertices[1] >= pntsForTri.size() ||
                t.vertices[2] >= pntsForTri.size())
                continue;

            std::vector<Vec3d> tri = {
                pntsForTri[t.vertices[0]],
                pntsForTri[t.vertices[1]],
                pntsForTri[t.vertices[2]]
            };

            // 如果原始是 CW（我们反转过），现在要再反转三角形以匹配原始法向
            if (needReverse) {
                std::reverse(tri.begin(), tri.end());
            }

            result.push_back(std::move(tri));
        }

        return result;
    }

    std::vector<std::vector<Vec3d>> GeomCalc::TriangulateWithConstraints(
        const std::vector<Vec3d>& outerBoundary,
        const std::vector<std::vector<Vec3d>>& intersectionPolylines,
        std::vector<Vec3d>& newIntPnts)
    {
        newIntPnts.clear();
        if (outerBoundary.size() < 3) return {};

        // === 1. 规范化外边界（移除重复首尾点）===
        std::vector<Vec3d> cdtOuter = outerBoundary;
        while (cdtOuter.size() > 2 && (cdtOuter.front() - cdtOuter.back()).LengthSq() <= g_epsilon * g_epsilon)
            cdtOuter.pop_back();
        if (cdtOuter.size() < 3) return {};

        // === 2. 计算投影平面 ===
        Vec3d origin = std::accumulate(cdtOuter.begin(), cdtOuter.end(), Vec3d{ 0,0,0 }) * (1.0 / cdtOuter.size());
        Vec3d normal = CompuateNormal(cdtOuter);
        if (normal.LengthSq() < g_epsilon * g_epsilon) return {};

        Vec3d x_axis = normal.Cross(std::abs(normal.x) < 0.9 ? Vec3d{ 1,0,0 } : Vec3d{ 0,1,0 }).Normalization();
        Vec3d y_axis = normal.Cross(x_axis);

        auto project = [&](const Vec3d& p) {
            Vec3d rel = p - origin;
            return CDT::V2d<double>{rel.Dot(x_axis), rel.Dot(y_axis)};
            };

        // === 3. 确保外边界为CCW（2D投影）===
        std::vector<CDT::V2d<double>> outer2D(cdtOuter.size());
        std::transform(cdtOuter.begin(), cdtOuter.end(), outer2D.begin(), project);

        double area = 0.0;
        for (size_t i = 0, n = outer2D.size(); i < n; ++i)
            area += outer2D[i].x * outer2D[(i + 1) % n].y - outer2D[(i + 1) % n].x * outer2D[i].y;

        if (area < 0) std::reverse(cdtOuter.begin(), cdtOuter.end());

        // === 4. 收集所有输入点（3D去重）并构建2D顶点 ===
        std::map<Vec3d, size_t, Vec3dCmp> pointToIndex;
        auto addPoint = [&](const Vec3d& p) {
            auto [it, inserted] = pointToIndex.emplace(p, pointToIndex.size());
            return it->second;
            };

        // 添加外边界点
        for (const auto& p : cdtOuter) addPoint(p);
        // 添加交线点
        for (const auto& poly : intersectionPolylines)
            for (const auto& p : poly) addPoint(p);

        // 构建CDT顶点（2D投影）
        std::vector<CDT::V2d<double>> cdtVertices(pointToIndex.size());
        for (const auto& [p, idx] : pointToIndex)
            cdtVertices[idx] = project(p);

        // === 5. 构建约束边（自动去重）===
        CDT::EdgeUSet edgeSet;
        auto addEdge = [&](size_t i, size_t j) {
            if (i != j) edgeSet.insert({ static_cast<CDT::VertInd>(i), static_cast<CDT::VertInd>(j) });
            };

        // 外边界（闭合）
        for (size_t i = 0, n = cdtOuter.size(); i < n; ++i)
            addEdge(pointToIndex[cdtOuter[i]], pointToIndex[cdtOuter[(i + 1) % n]]);

        // 交线（开链）
        for (const auto& poly : intersectionPolylines)
            for (size_t i = 0, n = poly.size() - 1; i < n; ++i)
                addEdge(pointToIndex[poly[i]], pointToIndex[poly[i + 1]]);

        std::vector<CDT::Edge> constraints(edgeSet.begin(), edgeSet.end());

        // === 6. 执行CDT（启用自动解决相交）===
        CDT::Triangulation<double> tri(
            CDT::VertexInsertionOrder::Auto,
            CDT::IntersectingConstraintEdges::TryResolve,
            g_epsilon
        );

        tri.insertVertices(cdtVertices);
        const size_t origVertCount = tri.vertices.size(); // CDT 2D去重后的原始顶点数
        tri.insertEdges(constraints);
        tri.eraseOuterTrianglesAndHoles();

        // === 7. 重建完整3D顶点（含交点）并提取新顶点 ===
        auto to3D = [&](const CDT::V2d<double>& v) {
            return origin + x_axis * v.x + y_axis * v.y;
            };

        std::vector<Vec3d> allPoints3D;
        allPoints3D.reserve(tri.vertices.size());
        std::transform(tri.vertices.begin(), tri.vertices.end(), std::back_inserter(allPoints3D), to3D);

        if (origVertCount < allPoints3D.size())
            newIntPnts.assign(allPoints3D.begin() + origVertCount, allPoints3D.end());

        // === 8. 提取三角形（恢复原始绕向）===
        std::vector<std::vector<Vec3d>> result;
        for (const auto& t : tri.triangles) {
            if (t.vertices[0] == CDT::noVertex || t.vertices[1] == CDT::noVertex || t.vertices[2] == CDT::noVertex)
                continue;

            std::vector<Vec3d> tri3D = { allPoints3D[t.vertices[0]], allPoints3D[t.vertices[1]], allPoints3D[t.vertices[2]] };
            if (area < 0) std::reverse(tri3D.begin(), tri3D.end()); // 与原始外边界绕向一致
            result.push_back(std::move(tri3D));
        }

        return result;
    }

    // 主函数
    bool GeomCalc::TriangulateWithHoles(
        const std::vector<Vec3d>& bnd,
        const std::vector<std::vector<Vec3d>>& holes,
        std::vector<std::vector<Vec3d>>& outTriangles)
    {
        const double eps = g_epsilon;
        const double epsSq = eps * eps;
        outTriangles.clear();

        if (bnd.size() < 3) return false;

        // === 1. 用 bnd 拟合平面 ===
        Vec3d origin{ 0, 0, 0 };
        for (const auto& p : bnd) {
            origin = origin + p;
        }
        origin = origin * (1.0 / static_cast<double>(bnd.size()));
        Vec3d refNormal = CompuateNormal(bnd);

        // === 2. 构建局部坐标系 ===
        Vec3d arbitrary = (std::abs(refNormal.x) < 0.9) ? Vec3d{ 1, 0, 0 } : Vec3d{ 0, 1, 0 };
        Vec3d x_axis = refNormal.Cross(arbitrary);
        if (x_axis.LengthSq() < epsSq) {
            arbitrary = Vec3d{ 0, 0, 1 };
            x_axis = refNormal.Cross(arbitrary);
        }
        x_axis = x_axis.Normalization();
        Vec3d y_axis = refNormal.Cross(x_axis);

        // === 3. 投影函数 ===
        auto to2D = [&](const Vec3d& p) -> PointD {
            Vec3d rel = p - origin;
            return PointD(rel.Dot(x_axis), rel.Dot(y_axis));
            };

        auto to3D = [&](const PointD& pt) -> Vec3d {
            return origin + x_axis * pt.x + y_axis * pt.y;
            };

        // === 4. 投影外轮廓并转为 CCW ===
        PathD outer = [&]() {
            PathD p;
            p.reserve(bnd.size());
            for (const auto& v : bnd) p.push_back(to2D(v));
            if (SignedArea(p) < 0) std::reverse(p.begin(), p.end());
            return p;
            }();

        // === 5. 投影洞，并转为 CW（相对于 outer 的 CCW）===
        PathsD innerPaths;
        for (const auto& hole : holes) {
            if (hole.size() < 3) continue;
            PathD h;
            h.reserve(hole.size());
            for (const auto& v : hole) h.push_back(to2D(v));
            // 洞应为 CW（即面积 < 0），若不是则反转
            if (SignedArea(h) > 0) std::reverse(h.begin(), h.end());
            innerPaths.push_back(h);
        }

        // === 6. 使用 Clipper2 执行布尔差集：outer \ union(holes) ===
        PathsD subjects = { outer };
        PathsD clips = innerPaths;

        // 先对 holes 做 Union（避免重叠洞导致问题）
        if (!clips.empty()) {
            clips = Union(clips, FillRule::EvenOdd);
        }

        PathsD diff;
        if (clips.empty()) {
            diff = subjects;
        }
        else {
            diff = Difference(subjects, clips, FillRule::EvenOdd);
        }

        if (diff.empty()) {
            return true; // 无剩余区域
        }

        // === 7. 提取所有环（外环+内环）用于 CDT ===
        // CDT 要求：顶点列表 + 约束边（包括外环和内环的边）
        std::vector<CDT::V2d<double>> cdtVertices;
        std::vector<CDT::Edge> cdtEdges;

        // 映射：2D点 -> 顶点索引（去重）
        struct  PDCompare {
            bool operator()(const PointD& pd1, const PointD& pd2) const {
                if (std::fabs(pd1.x - pd2.x) > g_epsilon)
                    return pd1.x < pd2.x;
                if (std::fabs(pd1.y - pd2.y) > g_epsilon)
                    return pd1.y < pd2.y;
                return false;
            }
        };
        std::map<PointD, std::size_t, PDCompare> pointToIndex;

        auto addPathToCDT = [&](const PathD& path) {
            if (path.size() < 3) return;
            std::vector<std::size_t> indices;
            for (const auto& pt : path) {
                auto it = pointToIndex.find(pt);
                if (it == pointToIndex.end()) {
                    std::size_t idx = cdtVertices.size();
                    cdtVertices.emplace_back(pt.x, pt.y);
                    pointToIndex[pt] = idx;
                    indices.push_back(idx);
                }
                else {
                    indices.push_back(it->second);
                }
            }
            // 添加约束边（闭合）
            for (size_t i = 0; i < indices.size(); ++i) {
                std::size_t j = (i + 1) % indices.size();
                cdtEdges.emplace_back(static_cast<CDT::VertInd>(indices[i]),
                    static_cast<CDT::VertInd>(indices[j]));
            }
            };

        for (const auto& poly : diff) {
            addPathToCDT(poly);
        }

        if (cdtVertices.empty()) return true;

        // === 8. 执行 CDT 三角化 ===
        CDT::Triangulation<double> cdt;
        cdt.insertVertices(cdtVertices);
        cdt.insertEdges(cdtEdges);
        cdt.eraseOuterTrianglesAndHoles();

        // === 9. 转换回 3D，并校正法向 ===
        outTriangles.clear();
        for (const auto& tri : cdt.triangles) {
            if (tri.vertices[0] == tri.vertices[1] ||
                tri.vertices[1] == tri.vertices[2] ||
                tri.vertices[0] == tri.vertices[2]) continue;

            Vec3d v0 = to3D(PointD(cdtVertices[tri.vertices[0]].x, cdtVertices[tri.vertices[0]].y));
            Vec3d v1 = to3D(PointD(cdtVertices[tri.vertices[1]].x, cdtVertices[tri.vertices[1]].y));
            Vec3d v2 = to3D(PointD(cdtVertices[tri.vertices[2]].x, cdtVertices[tri.vertices[2]].y));

            // 检查法向
            Vec3d n = (v1 - v0).Cross(v2 - v0);
            if (n.Dot(refNormal) < 0) {
                // 反转顺序
                std::swap(v1, v2);
            }

            outTriangles.push_back({ v0, v1, v2 });
        }

        return true;
    }


    // 辅助：提取三角形索引（处理 1-based + 0 分隔）
    static std::vector<std::array<unsigned int, 3>> extractTriangles(const TriMesh& mesh) {
        std::vector<std::array<unsigned int, 3>> tris;
        const auto& idx = mesh.indices;
        const size_t n = mesh.points.size();

        for (size_t i = 0; i < idx.size(); ) {
            if (idx[i] == 0) {
                ++i;
                continue;
            }
            if (i + 2 >= idx.size()) break;

            int i0 = idx[i] - 1;
            int i1 = idx[i + 1] - 1;
            int i2 = idx[i + 2] - 1;

            if (i0 >= 0 && i1 >= 0 && i2 >= 0 &&
                static_cast<size_t>(i0) < n &&
                static_cast<size_t>(i1) < n &&
                static_cast<size_t>(i2) < n) {
                tris.push_back({ static_cast<unsigned int>(i0),
                                static_cast<unsigned int>(i1),
                                static_cast<unsigned int>(i2) });
            }
            i += 4; // 跳过 3 个顶点 + 1 个分隔符
        }
        return tris;
    }

    // 判断是否为闭合实体
    bool GeomCalc::IsClosedSolid(const TriMesh& mesh) {
        if (mesh.points.empty() || mesh.indices.empty()) {
            return false;
        }

        auto triangles = extractTriangles(mesh);
        if (triangles.empty()) {
            return false;
        }

        // 使用 unordered_map 存储边 -> 出现次数
        // 边表示为 (min, max) 的 pair，确保方向无关
        struct EdgeHash {
            std::size_t operator()(const std::pair<unsigned int, unsigned int>& e) const {
                return (static_cast<std::size_t>(e.first) << 16) ^ e.second;
            }
        };

        std::unordered_map<std::pair<unsigned int, unsigned int>, int, EdgeHash> edgeCount;

        for (const auto& tri : triangles) {
            unsigned int a = tri[0], b = tri[1], c = tri[2];

            // 三条边：(a,b), (b,c), (c,a)
            std::pair<unsigned int, unsigned int> e1 = { std::min(a, b), std::max(a, b) };
            std::pair<unsigned int, unsigned int> e2 = { std::min(b, c), std::max(b, c) };
            std::pair<unsigned int, unsigned int> e3 = { std::min(c, a), std::max(c, a) };

            edgeCount[e1]++;
            edgeCount[e2]++;
            edgeCount[e3]++;
        }

        // 检查：每条边必须出现恰好 2 次
        for (const auto& kv : edgeCount) {
            if (kv.second != 2) {
                return false;
            }
        }

        return true;
    }