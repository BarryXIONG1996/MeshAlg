#include "BooleanOpHelper.h"

BooleanOpHelper::BooleanOpHelper(TopoTriMesh& objMesh, TopoTriMesh& subMesh, int opType)
    : m_objMesh(objMesh), m_subMesh(subMesh), m_opType(opType)
{
}

bool BooleanOpHelper::Execute(TopoTriMesh& res)
{
    // 根据m_opType执行相应的布尔操作
    switch (m_opType)
    {
    case 0:
        // 执行并集操作
        break;
    case 1:
        // 执行交集操作
        break;
    case 2:
        // 执行差集操作
        break;
    default:
        // 处理未知操作类型
        return false;
    }

    return true;
}