#pragma once
#include <vector>
#include <map>
#include <Geometry.h>

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

