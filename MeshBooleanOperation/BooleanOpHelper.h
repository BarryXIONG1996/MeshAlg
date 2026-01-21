#pragma once
#include "TopoTriMesh.h"
#include <set>
#include <queue>

/*Note: 注意修改工具集的面片法向*/

// 布尔运算后，未使用的拓扑会被释放
class BooleanOpHelper
{
 public:
    BooleanOpHelper(TopoTriMesh& objMesh, TopoTriMesh& subMesh, TopoTriMesh& coPlanes, int opType);

    bool Execute(TopoTriMesh& res);

private:
    void BFSExtractRegion(const TopoTriMesh& mesh, TopoTriMesh& Mout, TopoTriMesh& Min); // 共享交点/交线拓扑
    void AddFace(Face* f, TopoTriMesh& M, std::set<Vertex*>& vv, std::set<Edge*>& ve, std::set<Face*>& vf);
    void AddEdge(Edge* e, TopoTriMesh& M, std::set<Vertex*>& vv, std::set<Edge*>& ve);
    void AddAdjacentEdges(Face* f, Edge* skipEdge, std::queue<Edge*>& es, std::set<Edge*>& visitedE);

    void CombineTopoTriMesh(TopoTriMesh& M, std::vector<TopoTriMesh>& Ms);
    void ReleaseMeshExceptBoundary(TopoTriMesh& M);

private:
    TopoTriMesh& m_objMesh;
    TopoTriMesh& m_subMesh;
    TopoTriMesh& m_coPlanes;
    int m_opType; // 0-UNION, 1-INTERSECTION, 2-DIFFERENCE
};

