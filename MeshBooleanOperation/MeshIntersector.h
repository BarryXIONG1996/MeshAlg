#pragma once
#include "TopoTriMesh.h"

class MeshIntersector
{
public:
    MeshIntersector(TopoTriMesh& m1, TopoTriMesh& m2);

    // 执行网格求交操作
    bool Execute(TopoTriMesh& coPlanes);

private:
    TopoTriMesh& m_mesh1;
    TopoTriMesh& m_mesh2;
};

