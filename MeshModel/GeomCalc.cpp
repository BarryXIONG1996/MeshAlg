#include "GeomCalc.h"
#include "poly2tri/poly2tri.h"
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

double GeomCalc::Point2PlaneDistatnce(Vec3d const& pnt, Vec3d const& o, Vec3d const& n)
{
    return std::fabs(n.Dot(Vec3d(pnt - o)));
}

bool GeomCalc::TriRegionSplit(std::vector<Vec3d> const& tri1, std::vector<Vec3d> const& tri2, std::vector<std::vector<Vec3d>>& tri1OutTri2, std::vector<std::vector<Vec3d>>& tri2OutTri1, std::vector<Vec3d>& tri1CollapseTri2)
{
    return false;
}

std::vector<std::vector<Vec3d>> GeomCalc::Triangulate(std::vector<Vec3d> const& inputPnts)
{
    if (inputPnts.empty()) return {};

    // ===== 1. 清理重复点（保留原始索引）=====
    std::vector<Vec3d> cleaned;
    std::vector<size_t> origIndices;
    cleaned.push_back(inputPnts[0]);
    origIndices.push_back(0);

    for (size_t i = 1; i < inputPnts.size(); ++i) {
        if ((inputPnts[i] - cleaned.back()).Length() > g_epsilon) {
            cleaned.push_back(inputPnts[i]);
            origIndices.push_back(i);
        }
    }

    if (cleaned.size() > 1 && (cleaned.front() - cleaned.back()).Length() < g_epsilon) {
        cleaned.pop_back();
        origIndices.pop_back();
    }
    if (cleaned.size() < 3) return {};

    // ===== 2. 移除共线点（保守）=====
    auto removeCollinear = [&](std::vector<Vec3d>& pts, std::vector<size_t>& idxs) {
        bool changed;
        do {
            changed = false;
            std::vector<Vec3d> newPts;
            std::vector<size_t> newIdxs;
            size_t n = pts.size();
            if (n < 3) break;

            for (size_t i = 0; i < n; ++i) {
                const auto& prev = pts[(i - 1 + n) % n];
                const auto& curr = pts[i];
                const auto& next = pts[(i + 1) % n];
                Vec3d e1 = curr - prev;
                Vec3d e2 = next - curr;
                if (e1.Cross(e2).Length() < g_epsilon) {
                    changed = true;
                    continue;
                }
                newPts.push_back(curr);
                newIdxs.push_back(idxs[i]);
            }
            if (changed && newPts.size() >= 3) {
                pts = std::move(newPts);
                idxs = std::move(newIdxs);
            }
            else {
                break;
            }
        } while (true);
        };

    std::vector<Vec3d> pnts = cleaned;
    std::vector<size_t> indices = origIndices;
    removeCollinear(pnts, indices);
    if (pnts.size() < 3) return {};

    // ===== 3. 计算法向量 =====
    Vec3d normal{ 0, 0, 0 };
    size_t n = pnts.size();
    for (size_t i = 0; i < n; ++i) {
        const auto& curr = pnts[i];
        const auto& next = pnts[(i + 1) % n];
        normal.x += (curr.y - next.y) * (curr.z + next.z);
        normal.y += (curr.z - next.z) * (curr.x + next.x);
        normal.z += (curr.x - next.x) * (curr.y + next.y);
    }
    if (normal.Length() < g_epsilon) return {}; // 退化

    // ===== 4. 投影到最佳 2D 平面 =====
    enum Plane { XY, XZ, YZ };
    Plane projPlane;
    double ax = std::fabs(normal.x), ay = std::fabs(normal.y), az = std::fabs(normal.z);
    if (ax >= ay && ax >= az) projPlane = YZ;
    else if (ay >= az) projPlane = XZ;
    else projPlane = XY;

    // 创建 2D 点（使用 unique_ptr 自动管理）
    std::vector<std::unique_ptr<p2t::Point>> pointStorage;
    std::vector<p2t::Point*> polyline;
    pointStorage.reserve(pnts.size());
    polyline.reserve(pnts.size());

    for (const auto& p : pnts) {
        double u, v;
        switch (projPlane) {
        case XY: u = p.x; v = p.y; break;
        case XZ: u = p.x; v = p.z; break;
        case YZ: u = p.y; v = p.z; break;
        }
        pointStorage.emplace_back(std::make_unique<p2t::Point>(u, v));
        polyline.push_back(pointStorage.back().get());
    }

    // ===== 5. 强制逆时针（使用 2D 点）=====
    double area = 0.0;
    for (size_t i = 0; i < polyline.size(); ++i) {
        size_t j = (i + 1) % polyline.size();
        area += polyline[i]->x * polyline[j]->y - polyline[j]->x * polyline[i]->y;
    }
    if (area < 0) {
        // 顺时针 → 反转
        std::reverse(polyline.begin(), polyline.end());
        std::reverse(indices.begin(), indices.end());
    }

    // ===== 6. 调用 Poly2Tri =====
    try {
        p2t::CDT cdt(polyline); // 正确传入 Point*
        cdt.Triangulate();

        const auto& triangles = cdt.GetTriangles();
        std::vector<std::vector<Vec3d>> result;
        result.reserve(triangles.size());

        // 构建 2D -> 3D 映射（通过索引）
        for (auto* tri : triangles) {
            std::vector<Vec3d> face;
            face.reserve(3);
            for (int i = 0; i < 3; ++i) {
                p2t::Point* pt = tri->GetPoint(i);
                // 找到该 Point* 在 polyline 中的位置 → 得到索引 → 得到原始 3D 点
                auto it = std::find(polyline.begin(), polyline.end(), pt);
                if (it != polyline.end()) {
                    size_t localIdx = std::distance(polyline.begin(), it);
                    size_t origInputIdx = indices[localIdx];
                    face.push_back(inputPnts[origInputIdx]);
                }
                // 理论上 always found，因为 CDT 只用边界点
            }
            if (face.size() == 3) {
                // 可选：检查是否退化
                Vec3d e1 = face[1] - face[0];
                Vec3d e2 = face[2] - face[0];
                if (e1.Cross(e2).Length() >= g_epsilon) {
                    result.push_back(std::move(face));
                }
            }
        }
        return result;
    }
    catch (...) {
        return {};
    }
}
