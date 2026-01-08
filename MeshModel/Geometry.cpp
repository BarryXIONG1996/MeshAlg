#include "Geometry.h"

const double g_epsilon = 1e-6;

Vec3d Vec3d::operator-(const Vec3d& v) {
    return Vec3d{ x - v.x, y - v.y, z - v.z };
}

Vec3d Vec3d::operator+(const Vec3d& v) {
    return Vec3d{ x + v.x, y + v.y, z + v.z };
}

Vec3d Vec3d::operator*(double s) const {
    return Vec3d{ x * s, y * s, z * s };
}

double Vec3d::Dot(const Vec3d& v) const {
    return x * v.x + y * v.y + z * v.z;
}

Vec3d Vec3d::Cross(const Vec3d& v) const {
    return Vec3d{
        y * v.z - z * v.y,
        z * v.x - x * v.z,
        x * v.y - y * v.x
    };
}

double Vec3d::Length() const {
    return sqrt(x * x + y * y + z * z);
}

bool Vec3dCmp::operator()(const Vec3d& a, const Vec3d& b) const {
    if (fabs(a.x - b.x) > g_epsilon) return a.x < b.x;
    if (fabs(a.y - b.y) > g_epsilon) return a.y < b.y;
    return a.z < b.z;
}

BndBox3d BndBox3d::Intersect(BndBox3d const& bnd) {
    BndBox3d res;
    res.lowerBnd.x = std::max(lowerBnd.x, bnd.lowerBnd.x);
    res.lowerBnd.y = std::max(lowerBnd.y, bnd.lowerBnd.y);
    res.lowerBnd.z = std::max(lowerBnd.z, bnd.lowerBnd.z);
    res.upperBnd.x = std::min(upperBnd.x, bnd.upperBnd.x);
    res.upperBnd.y = std::min(upperBnd.y, bnd.upperBnd.y);
    res.upperBnd.z = std::min(upperBnd.z, bnd.upperBnd.z);
    return res;
}

bool BndBox3d::IsOut(BndBox3d bnd) {
    // 判断两个包围盒是否不相交(精度g_epsilon)
    return (upperBnd.x < bnd.lowerBnd.x - g_epsilon || lowerBnd.x > bnd.upperBnd.x + g_epsilon ||
        upperBnd.y < bnd.lowerBnd.y - g_epsilon || lowerBnd.y > bnd.upperBnd.y + g_epsilon ||
        upperBnd.z < bnd.lowerBnd.z - g_epsilon || lowerBnd.z > bnd.upperBnd.z + g_epsilon);
}

void BndBox3d::Add(const Vec3d& p) {
    if (p.x < lowerBnd.x) lowerBnd.x = p.x;
    if (p.y < lowerBnd.y) lowerBnd.y = p.y;
    if (p.z < lowerBnd.z) lowerBnd.z = p.z;
    if (p.x > upperBnd.x) upperBnd.x = p.x;
    if (p.y > upperBnd.y) upperBnd.y = p.y;
    if (p.z > upperBnd.z) upperBnd.z = p.z;
}
