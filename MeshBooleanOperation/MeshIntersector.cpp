#include "MeshIntersector.h"
#include <memory>
#include <cmath>

// 定义 M_PI 和 M_E，防止未声明错误
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#include "RTree.h"

class SpatialAccelerator {
public:
    virtual void Build(std::vector<Face*>& faces) = 0;
    virtual std::vector<Face*> Query(const BndBox3d& box) = 0;
};

class RTreeAccelerator : public SpatialAccelerator {
    virtual void Build(std::vector<Face*>& faces) override 
    {
        // 构建R树
        for (Face* face : faces) {
            // 将每个面片的包围盒插入到R树中
            m_rtree.Insert((double*)&face->bbox.lowerBnd, (double*)&face->bbox.upperBnd, face);
        }
    }

    virtual std::vector<Face*> Query(const BndBox3d& box) override 
    {
        // 查询R树
        std::vector<Face*> result;
        m_rtree.Search((double*)&box.lowerBnd, (double*)&box.upperBnd,
            [&result](Face* face)->bool {
                result.push_back(face);
                return true; // 继续搜索
            });
        return result;
    }

    RTree<Face*, double, 3> m_rtree;
};

MeshIntersector::MeshIntersector(TopoTriMesh& m1, TopoTriMesh& m2)
    : m_mesh1(m1), m_mesh2(m2)
{

}

bool MeshIntersector::Execute(TopoTriMesh& coPlanes)
{
    // 1、构建空间搜索树
    // 计算工具集和目标集的包围盒
    BndBox3d boundingBoxO = m_mesh1.GetBndBox();
    BndBox3d boundingBoxT = m_mesh2.GetBndBox();
    // 计算包围盒的交集
    BndBox3d intBox = boundingBoxO.Intersect(boundingBoxT);
    // 遍历工具集中的所有三角面片
    std::vector<Face*>const& facesT = m_mesh2.fs;
    std::vector<Face*> intersectingFaces;
    for (Face* face : facesT) {
        // 检查三角面片的包围盒是否与intBox相交
        if (!face->bbox.IsOut(intBox)) {
            intersectingFaces.push_back(face);
        }
    }
    // 将相交的面片添加到空间搜索树
    std::shared_ptr<SpatialAccelerator> accelerator = std::make_shared<RTreeAccelerator>();
    accelerator->Build(intersectingFaces);

    // 2、相交测试
    std::vector<Face*>const& facesO = m_mesh1.fs;
    for (Face* face : facesO)
    {
        // 查询空间搜索树，获取可能相交的面片
        std::vector<Face*> candidateFaces = accelerator->Query(face->bbox);
        for (Face* candidateFace : candidateFaces) {
            // 进行精确的三角形相交计算
            // 如果相交，计算交点并更新网格数据结构

        }
    }
    
    return true;
}
