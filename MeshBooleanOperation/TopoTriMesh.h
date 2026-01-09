#pragma once
#include <vector>
#include <map>
#include <Geometry.h>

struct Edge;
struct Vertex {
    Vec3d pnt;
    Edge* e;        // 任一边
    int posTag;     // 0=default, 1=in, 2=out, 3=on
    std::vector<Edge*> GetAdjacentEdges();
};

struct Face;
struct Edge {
    Vertex* v1, * v2;
    Face* lF, * rF;  // 左右面
    Edge* lPE, * lSE, * rPE, * rSE; // 翼边结构
    bool isInner;

    std::vector<Vec3d> getPnts(bool left = true);
};

struct Face {
    Edge* e;        // 任一边
    BndBox3d bbox;

    std::vector<Edge*> getEdges();
    std::vector<Vec3d> getPnts();
};

struct TopoTriMesh{
    std::vector<Vertex*> vs;
    std::vector<Edge*> es;
    std::vector<Face*> fs;
    std::map<Vec3d, Vertex*, Vec3dCmp> p2V; // 几何去重（带epsilon）

    BndBox3d GetBndBox();
    void AddFace2TopoTriMesh(std::vector<Vec3d>const& pnts);
    void ReleaseMem();
};

