#pragma once
#include "TopoTriMesh.h"

class BooleanOpHelper
{
 public:
    BooleanOpHelper(TopoTriMesh& objMesh, TopoTriMesh& subMesh, int opType);

    bool Execute(TopoTriMesh& res);

private:
    TopoTriMesh& m_objMesh;
    TopoTriMesh& m_subMesh;
    int m_opType; // 0-UNION, 1-INTERSECTION, 2-DIFFERENCE
};

