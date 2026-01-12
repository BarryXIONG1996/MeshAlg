#pragma once
#include "TopoTriMesh.h"
#include <functional>

class MeshIntersector
{
public:
    MeshIntersector(TopoTriMesh& m1, TopoTriMesh& m2);
    // 执行网格求交操作
    bool Execute(TopoTriMesh& coPlanes);

private:
    std::vector<Vec3d> EdgeEdgeInt(Edge* e1, Edge* e2); // 返回值数组大小：0-没有交点,1-一个交点,2-重叠
    void FaceFaceInt(Face* f1, Face* f2, std::function<void(Face*, Face*)> intersectFunc, bool coPlanar);
    void CoPlanarFaceInt(Face* f1, Face* f2);
    void NonCoPlanarFaceInt(Face* f1, Face* f2);

private:
    TopoTriMesh& m_mesh1;
    TopoTriMesh& m_mesh2;
    TopoTriMesh m_coPlanes;

    std::map<std::pair<Edge*, Face*>, std::set<int>> m_ef2Int; // 边面对:交点
    std::vector<Face*, std::set<int>> m_face2IntSegs; // 面:交线
    std::vector<Vec3d> m_intersectPnts; // 交点
    std::vector<double> m_weights; // 交点权重
    std::vector<std::pair<int, int>> m_intSegs; // 交线
};