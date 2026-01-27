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

extern void ShowTriMesh(const TriMesh& mesh);

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

#ifdef _DEBUG
    //TriMesh t, o, co;
    //m_obj->ToMesh(o);
    //m_sub->ToMesh(t);
    //coPlanes.ToMesh(co);
    //ShowTriMesh(o);
    //ShowTriMesh(t);
    //ShowTriMesh(co);
#endif

    // 2.布尔运算
    BooleanOpHelper boolHelper(*m_obj, *m_sub, coPlanes, m_opType);
    if (!boolHelper.Execute(res))
        return false;

#ifdef _DEBUG
    TriMesh re;
    res.ToMesh(re);
    ShowTriMesh(re);
#endif

    return true;
}
