#include "MeshIntersector.h"
#include <GeomCalc.h>
#include <memory>
#include <cmath>
#include <set>

// 定义 M_PI 和 M_E，防止未声明错误
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#include "RTree.h"

extern const double g_epsilon;

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
        // 过滤掉与包围盒不相交的三角面片
        if (face->bbox.IsOut(intBox))
            continue;
        // 查询空间搜索树，获取可能相交的面片
        std::vector<Face*> candidateFaces = accelerator->Query(face->bbox);
        for (Face* candidateFace : candidateFaces) {
            // 3、面片之间的相交计算
            FaceFaceInt(face, candidateFace);
        }
    }
    
    coPlanes = m_coPlanes;

    return true;
}

std::vector<Vec3d> MeshIntersector::EdgeEdgeInt(Edge* e1, Edge* e2)
{
    return std::vector<Vec3d>();
}

void MeshIntersector::FaceFaceInt(Face* f1, Face* f2)
{
    if (f1->bbox.IsOut(f2->bbox))
        return;

    std::vector<Vec3d> ps1 = f1->getPnts(), ps2 = f2->getPnts();
    if (ps1.size() < 3 || ps2.size() < 3)
        return;

    Vec3d n1 = GeomCalc::CompuateNormal(ps1), n2 = GeomCalc::CompuateNormal(ps2);

    if (n1.Length() < g_epsilon || n2.Length() < g_epsilon)
        return;

    if (n1.Parallel(n2)) // 面平行
    {
        if (GeomCalc::Point2PlaneDistatnce(ps2.at(0), ps1.at(0), n1) > g_epsilon
            || GeomCalc::Point2PlaneDistatnce(ps1.at(0), ps2.at(0), n2) > g_epsilon)
            return; // 不共面

        // 共面
        CoPlanarFaceInt(f1, f2);

        return;
    }

    // 面不平行
    NonCoPlanarFaceInt(f1, f2);

    return;
}

void MeshIntersector::CoPlanarFaceInt(Face* f1, Face* f2)
{
    
}

void MeshIntersector::NonCoPlanarFaceInt(Face* f1, Face* f2)
{
}
