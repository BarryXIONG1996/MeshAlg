#include "BooleanOperation.h"
#include "MeshIntersector.h"
#include "BooleanOpHelper.h"

const double g_epsilon = 1e-9;

BooleanOperation::BooleanOperation(
    TopoTriMesh* operandA,
    TopoTriMesh* operandB,
    BooleanType op
) : m_obj(operandA), m_sub(operandB), m_opType(op)
{
}

 bool BooleanOperation::Execute(TopoTriMesh& res)
{
     if (!m_obj || !m_sub)
         return false;

    TopoTriMesh result;
    
    // 1. 求交
    MeshIntersector intersector(*m_obj, *m_sub);
    TopoTriMesh coPlanes;
    if (!intersector.Execute(coPlanes))
        return false;

    // 2.布尔运算
    BooleanOpHelper boolHelper(*m_obj, *m_sub, coPlanes, m_opType);
    if (!boolHelper.Execute(res))
        return false;

    return true;
}
