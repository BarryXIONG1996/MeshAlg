#pragma once
#include <vector>
#include <map>
#include <set>
#include <Geometry.h>

struct Edge;
struct Face;
struct Vertex {
    Vec3d pnt;
    Edge* e;        // 任一边
    int posTag;     // 0=default, 1=in, 2=out, 3=on
    std::vector<Edge*> GetAdjacentEdges();
    std::vector<Face*> GetAdjacentFaces();
};

struct Face;
struct Edge {
    Vertex* v1, * v2;
    Face* lF, * rF;  // (前进方向的)左右面
    Edge* lPE, * lSE, * rPE, * rSE; // (法线方向的)前后边
    bool isInner;

    std::vector<Vec3d> getPnts(bool left = true);
};

struct TopoTriMesh;
struct Face {
    Edge* es[3];
    bool reverse[3];
    BndBox3d bbox;
    TopoTriMesh* topo;

    std::vector<Vertex*> getVertices();
    std::vector<Edge*> getEdges();
    std::vector<bool> getEdgeDir();
    std::vector<Vec3d> getPnts();
};

struct TopoTriMesh{
    std::vector<Vertex*> vs;
    std::vector<Edge*> es;
    std::vector<Face*> fs;
    std::map<Vec3d, Vertex*, Vec3dCmp> p2V; // 几何去重（带epsilon）

    BndBox3d GetBndBox();
    Face* AddFace2TopoTriMesh(std::vector<Vec3d>const& pnts);
    void RemoveFace(Face* f);
    void RemoveEdge(Edge* e); // 移除边之前，要先移除两侧面
    void ReleaseMem();
    void Build(TriMesh const& triMesh);
    void ToMesh(TriMesh& mesh);
};

