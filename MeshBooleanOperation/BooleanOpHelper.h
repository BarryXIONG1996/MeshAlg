#pragma once
#include "TopoTriMesh.h"
#include <set>
#include <queue>

class BooleanOpHelper
{
 public:
    BooleanOpHelper(TopoTriMesh& objMesh, TopoTriMesh& subMesh, TopoTriMesh& coPlanes, int opType);

    bool Execute(TopoTriMesh& res);

private:
    void BFSExtractRegion(const TopoTriMesh& mesh, TopoTriMesh& Mout, TopoTriMesh& Min);
    void AddFace(Face* f, TopoTriMesh& M, std::set<Vertex*>& vv, std::set<Edge*>& ve, std::set<Face*>& vf);
    void AddEdge(Edge* e, TopoTriMesh& M, std::set<Vertex*>& vv, std::set<Edge*>& ve);
    void AddAdjacentEdges(Face* f, Edge* skipEdge, std::queue<Edge*>& es, std::set<Edge*>& visitedE);

    TopoTriMesh CombineMeshes(const TopoTriMesh& mesh1, const TopoTriMesh& mesh2, const TopoTriMesh& mesh3);

private:
    TopoTriMesh& m_objMesh;
    TopoTriMesh& m_subMesh;
    TopoTriMesh& m_coPlanes;
    int m_opType; // 0-UNION, 1-INTERSECTION, 2-DIFFERENCE
};

