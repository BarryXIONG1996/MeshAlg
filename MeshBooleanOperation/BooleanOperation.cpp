#include "BooleanOperation.h"
#include "MeshIntersector.h"
#include "BooleanOpHelper.h"

BooleanOperation::BooleanOperation(
    const TopoTriMesh& operandA,
    const TopoTriMesh& operandB,
    BooleanType op
) : m_obj(operandA), m_sub(operandB), m_opType(op)
{
}

 bool BooleanOperation::Execute(TopoTriMesh& res)
{
    TopoTriMesh result;
    
    // 1. 求交
    MeshIntersector intersector(m_obj, m_sub);
    if (!intersector.Execute())
        return false;

    // 2.布尔运算
    BooleanOpHelper boolHelper(m_obj, m_sub, m_opType);
    if (!boolHelper.Execute(res))
        return false;

    return true;
}
