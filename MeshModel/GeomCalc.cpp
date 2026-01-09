#include "GeomCalc.h"
#include <cmath>

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
