#include "Geometry.h"

const double g_epsilon = 1e-9;

Vec2d Vec2d::operator+(const Vec2d& v) const {
    return { x + v.x, y + v.y };
}

Vec2d Vec2d::operator-(const Vec2d& v) const {
    return { x - v.x, y - v.y };
}

Vec2d Vec2d::operator*(double s) const {
    return { x * s, y * s };
}

Vec2d Vec2d::operator/(double s) const {
    return { x / s, y / s };
}

// 点积
double Vec2d::Dot(const Vec2d& v) const {
    return x * v.x + y * v.y;
}

// 2D 叉积（标量）
double Vec2d::Cross(const Vec2d& v) const {
    return x * v.y - y * v.x;
}

// 长度
double Vec2d::Length() const {
    return sqrt(x * x + y * y);
}

// 长度平方
double Vec2d::LengthSq() const {
    return x * x + y * y;
}

// 归一化
Vec2d Vec2d::Normalization() const {
    double len = Length();
    if (len  < g_epsilon) {
        return { 0.0, 0.0 };
    }
    return { x / len, y / len };
}

// 平行判断
bool Vec2d::Parallel(const Vec2d& v) const {
    return std::abs(Cross(v)) < g_epsilon;
}

// 相等判断
bool Vec2d::Equal(const Vec2d& v) const {
    return operator-(v).Length() < g_epsilon;
}

Vec3d Vec3d::operator-(const Vec3d& v) const {
    return Vec3d{ x - v.x, y - v.y, z - v.z };
}

Vec3d Vec3d::operator+(const Vec3d& v) const {
    return Vec3d{ x + v.x, y + v.y, z + v.z };
}

Vec3d Vec3d::operator*(double s) const {
    return Vec3d{ x * s, y * s, z * s };
}

Vec3d Vec3d::operator/(double s) const {
    return { x / s, y / s, z / s };
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

double Vec3d::LengthSq() const
{
    return x * x + y * y + z * z;
}

Vec3d Vec3d::Normalization() const
{
    double len = Length();

    if (len < g_epsilon)
        return { 0,0,0 };

    return { x/len, y/len, z/len };
}

bool Vec3d::Parallel(const Vec3d& v) const
{
    Vec3d crossVec = Cross(v);
    return crossVec.Length() < g_epsilon;
}

bool Vec3d::Equal(const Vec3d& v) const
{
    return operator-(v).Length() < g_epsilon;
}

bool Vec3dCmp::operator()(const Vec3d& a, const Vec3d& b) const {
    if (fabs(a.x - b.x) > g_epsilon) return a.x < b.x;
    if (fabs(a.y - b.y) > g_epsilon) return a.y < b.y;
    return a.z < b.z - g_epsilon;
}

BndBox3d::BndBox3d()
{
    lowerBnd = { 
        std::numeric_limits<double>::max(), 
        std::numeric_limits<double>::max(), 
        std::numeric_limits<double>::max() 
    };

    upperBnd = {
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest()
    };
}

BndBox3d BndBox3d::Intersect(BndBox3d const& bnd) const {
    BndBox3d res;
    res.lowerBnd.x = std::max(lowerBnd.x, bnd.lowerBnd.x);
    res.lowerBnd.y = std::max(lowerBnd.y, bnd.lowerBnd.y);
    res.lowerBnd.z = std::max(lowerBnd.z, bnd.lowerBnd.z);
    res.upperBnd.x = std::min(upperBnd.x, bnd.upperBnd.x);
    res.upperBnd.y = std::min(upperBnd.y, bnd.upperBnd.y);
    res.upperBnd.z = std::min(upperBnd.z, bnd.upperBnd.z);
    return res;
}

bool BndBox3d::IsOut(BndBox3d bnd) const {
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

Vec3d BndBox3d::Center() const {
    return {
        (lowerBnd.x + upperBnd.x) / 2,
        (lowerBnd.y + upperBnd.y) / 2,
        (lowerBnd.z + upperBnd.z) / 2,
    };
}
