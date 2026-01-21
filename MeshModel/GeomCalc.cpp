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
bool GeomCalc::IsLeft(const Vec3d& a, const Vec3d& b, const Vec3d& p)
{
    Vec3d ab = b - a;
    Vec3d ap = p - a;

    // 计算叉积 (用于找最大投影面)
    Vec3d cross = ab.Cross(ap);

    // 找绝对值最大的分量（决定投影到哪个坐标平面）
    double ax = std::abs(cross.x);
    double ay = std::abs(cross.y);
    double az = std::abs(cross.z);

    double area2;
    if (ax >= ay && ax >= az) {
        // 投影到 YZ 平面（忽略 x）
        area2 = ab.y * ap.z - ab.z * ap.y;
    }
    else if (ay >= az) {
        // 投影到 XZ 平面（忽略 y）
        area2 = ab.z * ap.x - ab.x * ap.z;
    }
    else {
        // 投影到 XY 平面（忽略 z）
        area2 = ab.x * ap.y - ab.y * ap.x;
    }

    return area2 > g_epsilon;
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



// 主函数：多边形求交（共面假设）
bool GeomCalc::PolyIntersect(
    const std::vector<Vec3d>& poly1,
    const std::vector<Vec3d>& poly2,
    std::vector<Vec3d>& outIntersection)
{
    // 辅助：计算多边形有向面积（2D）
    auto SignedArea = [](const PathD& p) {
        if (p.size() < 3) return 0.0;
        double area = 0.0;
        size_t n = p.size();
        for (size_t i = 0; i < n; ++i) {
            size_t j = (i + 1) % n;
            area += p[i].x * p[j].y - p[j].x * p[i].y;
        }
        return area * 0.5;
    };

    const double eps = g_epsilon;
    const double epsSq = eps * eps;
    outIntersection.clear();

    if (poly1.size() < 3 || poly2.size() < 3) {
        return true; // 无效多边形，无交集
    }

    // === 1. 用 poly1 计算参考法向和平面 ===
    Vec3d origin = poly1[0];
    Vec3d normalSum{ 0, 0, 0 };
    bool hasValidNormal = false;

    // 遍历 poly1 的连续三角形扇区，累加法向
    for (size_t i = 1; i + 1 < poly1.size(); ++i) {
        Vec3d e1 = poly1[i] - poly1[0];
        Vec3d e2 = poly1[i + 1] - poly1[0];
        Vec3d n = e1.Cross(e2);
        if (n.LengthSq() > epsSq) {
            normalSum = normalSum + n;
            hasValidNormal = true;
        }
    }

    if (!hasValidNormal || normalSum.LengthSq() < epsSq) {
        // poly1 退化（共线或点）
        return true;
    }

    Vec3d refNormal = normalSum.Normalization();

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
        const std::vector<std::vector<Vec3d>>& intersectionPolylines)
    {
        if (outerBoundary.size() < 3) {
            return {};
        }

        // --- Step 1: 计算平面参数 ---
        Vec3d origin{ 0, 0, 0 };
        for (const auto& p : outerBoundary) {
            origin = origin + p;
        }
        origin = origin * (1.0 / static_cast<double>(outerBoundary.size()));

        Vec3d normal = CompuateNormal(outerBoundary);
        if (normal.LengthSq() < g_epsilon * g_epsilon) {
            return {}; // 退化
        }

        Vec3d arbitrary = (std::abs(normal.x) < 0.9) ? Vec3d{ 1, 0, 0 } : Vec3d{ 0, 1, 0 };
        Vec3d x_axis = normal.Cross(arbitrary).Normalization();
        Vec3d y_axis = normal.Cross(x_axis);

        // 投影函数
        auto project = [&](const Vec3d& p) -> CDT::V2d<double> {
            Vec3d rel = p - origin;
            double u = rel.Dot(x_axis);
            double v = rel.Dot(y_axis);
            return { u, v };
            };

        // --- Step 2: 投影外边界到 2D，判断绕向 ---
        std::vector<CDT::V2d<double>> outer2D;
        outer2D.reserve(outerBoundary.size());
        for (const auto& p : outerBoundary) {
            outer2D.push_back(project(p));
        }

        // 计算有符号面积（标准公式）
        double signedArea = 0.0;
        size_t N = outer2D.size();
        for (size_t i = 0; i < N; ++i) {
            size_t j = (i + 1) % N;
            signedArea += outer2D[i].x * outer2D[j].y - outer2D[j].x * outer2D[i].y;
        }

        bool needReverseOuter = (signedArea < 0); // <0 表示 CW（在右手系中）

        // 决定用于 CDT 的外边界（必须是 CCW）
        std::vector<Vec3d> cdtOuter = outerBoundary;
        if (needReverseOuter) {
            std::reverse(cdtOuter.begin(), cdtOuter.end());
        }

        // --- Step 3: 收集所有唯一点（使用原始几何，含反转后的外边界）---
        std::map<Vec3d, size_t, Vec3dCmp> pointToIndex;
        std::vector<Vec3d> uniquePoints;

        auto addPoint = [&](const Vec3d& p) -> size_t {
            auto it = pointToIndex.find(p);
            if (it != pointToIndex.end()) return it->second;
            size_t idx = uniquePoints.size();
            uniquePoints.push_back(p);
            pointToIndex[p] = idx;
            return idx;
            };

        // 添加外边界（可能是反转后的）
        for (const auto& p : cdtOuter) {
            addPoint(p);
        }

        // 添加交线（保持原始顺序）
        for (const auto& poly : intersectionPolylines) {
            for (const auto& p : poly) {
                addPoint(p);
            }
        }

        // --- Step 4: 构建 CDT 顶点 ---
        std::vector<CDT::V2d<double>> cdtVertices;
        cdtVertices.reserve(uniquePoints.size());
        for (const auto& p : uniquePoints) {
            cdtVertices.push_back(project(p));
        }

        // --- Step 5: 构建约束边 ---
        std::vector<CDT::Edge> constraints;

        // 外边界（现在是 CCW）
        for (size_t i = 0; i < cdtOuter.size(); ++i) {
            size_t j = (i + 1) % cdtOuter.size();
            CDT::VertInd vi = static_cast<CDT::VertInd>(pointToIndex.at(cdtOuter[i]));
            CDT::VertInd vj = static_cast<CDT::VertInd>(pointToIndex.at(cdtOuter[j]));
            constraints.emplace_back(vi, vj);
        }

        // 交线（开链，原始顺序）
        for (const auto& poly : intersectionPolylines) {
            for (size_t i = 0; i + 1 < poly.size(); ++i) {
                CDT::VertInd vi = static_cast<CDT::VertInd>(pointToIndex.at(poly[i]));
                CDT::VertInd vj = static_cast<CDT::VertInd>(pointToIndex.at(poly[i + 1]));
                constraints.emplace_back(vi, vj);
            }
        }

        // --- Step 6: 执行 CDT ---
        CDT::Triangulation<double> tri;
        tri.insertVertices(cdtVertices);
        tri.insertEdges(constraints);
        tri.eraseOuterTrianglesAndHoles();

        // --- Step 7: 提取结果，并恢复原始绕向 ---
        std::vector<std::vector<Vec3d>> result;
        for (const auto& t : tri.triangles) {
            if (t.vertices[0] == CDT::noVertex ||
                t.vertices[1] == CDT::noVertex ||
                t.vertices[2] == CDT::noVertex) {
                continue;
            }

            std::vector<Vec3d> triangle = {
                uniquePoints[t.vertices[0]],
                uniquePoints[t.vertices[1]],
                uniquePoints[t.vertices[2]]
            };

            // 如果原始外边界是 CW（我们反转过），则三角形也要反转以匹配原始法向
            if (needReverseOuter) {
                std::reverse(triangle.begin(), triangle.end());
            }

            result.push_back(std::move(triangle));
        }

        return result;
    }