#pragma once

#ifdef MESHMODEL_EXPORTS
#define MESHMODELDLL __declspec(dllexport)
#else
#define MESHMODELDLL __declspec(dllimport)
#endif

#include <vector>
#include <map>

struct MESHMODELDLL Vec3d
{
    double x, y, z;
    Vec3d operator-(const Vec3d& v);
    Vec3d operator+(const Vec3d& v);
    Vec3d operator*(double s) const;
    double Dot(const Vec3d& v) const;
    Vec3d Cross(const Vec3d& v) const;
    double Length() const;
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

