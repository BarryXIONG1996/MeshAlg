#pragma once
#include "TopoTriMesh.h"

class MeshIntersector
{
public:
    MeshIntersector(TopoTriMesh& m1, TopoTriMesh& m2);
    // 执行网格求交操作
    bool Execute(TopoTriMesh& coPlanes);

private:
    std::vector<Vec3d> EdgeEdgeInt(Edge* e1, Edge* e2); // 返回值数组大小：0-没有交点,1-一个交点,2-重叠
    void FaceFaceInt(Face* f1, Face* f2);
    void CoPlanarFaceInt(Face* f1, Face* f2);
    void NonCoPlanarFaceInt(Face* f1, Face* f2);

private:
    TopoTriMesh& m_mesh1;
    TopoTriMesh& m_mesh2;
    TopoTriMesh m_coPlanes;
};

