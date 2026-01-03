#pragma once
#include <vector>
#include <map>

struct Vec3d 
{ 
    double x, y, z;
    Vec3d operator-(const Vec3d& v);
    Vec3d operator+(const Vec3d& v);
    Vec3d operator*(double s) const;
    double Dot(const Vec3d& v) const;
    Vec3d Cross(const Vec3d& v) const;
    double Length() const;
};

struct Vec3dCmp {
    bool operator()(const Vec3d& a, const Vec3d& b) const;
};

struct BndBox3d 
{ 
    Vec3d lowerBnd, upperBnd;
    BndBox3d Intersect(BndBox3d const& bnd);
    bool IsOut(BndBox3d bnd);
    void Add(const Vec3d& p);
};

struct Edge;
struct Vertex {
    Vec3d pnt;
    Edge* e;        // 任一边
    int posTag;     // 0=default, 1=in, 2=out, 3=on
};

struct Face;
struct Edge {
    Vertex* v1, * v2;
    Face* lF, * rF;  // 左右面
    Edge* lPE, * lSE, * rPE, * rSE; // 翼边结构
    bool isInner;
};

struct Face {
    Edge* e;        // 任一边
    BndBox3d bbox;
};

struct TopoTriMesh{
    std::vector<Vertex*> vs;
    std::vector<Edge*> es;
    std::vector<Face*> fs;
    std::map<Vec3d, Vertex*, Vec3dCmp> p2V; // 几何去重（带epsilon）

    BndBox3d GetBndBox();
};

