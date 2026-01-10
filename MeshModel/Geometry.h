#pragma once

#ifdef MESHMODEL_EXPORTS
#define MESHMODELDLL __declspec(dllexport)
#else
#define MESHMODELDLL __declspec(dllimport)
#endif

#include <vector>
#include <map>

extern const double g_epsilon;

struct MESHMODELDLL Vec2d
{
    double x, y;
    Vec2d operator+(const Vec2d& v) const;
    Vec2d operator-(const Vec2d& v) const;
    Vec2d operator*(double s) const;
    Vec2d operator/(double s) const;
    double Dot(const Vec2d& v) const;
    double Cross(const Vec2d& v) const;
    double Length() const;
    double LengthSq() const;
    Vec2d Normalization() const;
    bool Parallel(const Vec2d& v) const;
    bool Equal(const Vec2d& v) const;
};

struct MESHMODELDLL Vec3d
{
    double x, y, z;
    Vec3d operator-(const Vec3d& v) const;
    Vec3d operator+(const Vec3d& v) const;
    Vec3d operator*(double s) const;
    Vec3d operator/(double s) const;
    double Dot(const Vec3d& v) const;
    Vec3d Cross(const Vec3d& v) const;
    double Length() const;
    double LengthSq() const;
    Vec3d Normalization() const;
    bool Parallel(const Vec3d& v) const;
    bool Equal(const Vec3d& v) const;
};

struct MESHMODELDLL Vec3dCmp {
    bool operator()(const Vec3d& a, const Vec3d& b) const;
};

struct MESHMODELDLL BndBox3d
{
    Vec3d lowerBnd, upperBnd;
    BndBox3d Intersect(BndBox3d const& bnd);
    bool IsOut(BndBox3d bnd);
    void Add(const Vec3d& p);
};

