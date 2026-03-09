#pragma once
#include "TopoTriMesh.h"

/*布尔算法类*/
class BooleanOperation
{
public:
    enum BooleanType { UNION = 0, INTERSECTION, DIFFERENCE };

    BooleanOperation(
        TopoTriMesh* operandA,   // 目标集
        TopoTriMesh* operandB,   // 工具集 
        BooleanType op           // 布尔操作类型：UNION / INTERSECTION / DIFFERENCE
    );

    bool Execute(TopoTriMesh& res);    // 执行布尔运算，返回结果集

private:
    TopoTriMesh* m_obj/*目标集*/, *m_sub/*工具集*/;
    BooleanType m_opType;   
};