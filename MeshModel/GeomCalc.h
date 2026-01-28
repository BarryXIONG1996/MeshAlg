#pragma once
#include "Geometry.h"

struct MESHMODELDLL GeomCalc
{
    // 计算面片法向
    static Vec3d CompuateNormal(std::vector<Vec3d>const& ps);

    // 判断点与线段的位置关系
    static bool IsLeft(const Vec3d& a, const Vec3d& b, const Vec3d& p);

    // 判断点是否在线段上
    static bool IsPointOnSegment(Vec3d const& p, Vec3d const& a, Vec3d const& b, double& paramP/*点p的参数*/);

    // 计算点到面的距离
    static double Point2PlaneDistatnce(Vec3d const& pnt, Vec3d const& o, Vec3d const& n);
    static double Point2PlaneSignDistance(Vec3d const& pnt, Vec3d const& o, Vec3d const& n);

    // 输出线段上的有序点
    static std::vector<Vec3d> OrderPointsOnSegment(
        const Vec3d& start,
        const Vec3d& end,
        const std::vector<Vec3d>& pts);

    // 计算面面交线
    static bool CalPlanePlaneIntersection(Vec3d const& o1, Vec3d const& dir1, Vec3d const& o2, Vec3d const& dir2, Vec3d& intO, Vec3d& intDir);

    // 计算直线和线段交点
    static bool CalLineSegmentIntersection(Vec3d const& o, Vec3d const& dir, Vec3d const& segS, Vec3d const& segE, double& lnP1/*直线参数*/, double& SegP2/*线段参数*/);

    // 直线和平面求交
    static bool linePlaneIntersect(
        const Vec3d& lineOrigin,
        const Vec3d& lineDir,
        const Vec3d& planeOrigin,
        const Vec3d& planeNormal,
        double& t);

    // 多边形求交
    static bool PolyIntersect(
        const std::vector<Vec3d>& poly1,
        const std::vector<Vec3d>& poly2,
        std::vector<Vec3d>& outIntersection);

    // 多边形的三角化（结果法向保持一致）
    static std::vector<std::vector<Vec3d>> Triangulate(std::vector<Vec3d> const& poly);

    // 包含限制边的三角化
    static std::vector<std::vector<Vec3d>> TriangulateWithConstraints(
        const std::vector<Vec3d>& outerBoundary,
        const std::vector<std::vector<Vec3d>>& intersectionPolylines);

    // 对内外环多边形进行三角化
    static bool TriangulateWithHoles(
        const std::vector<Vec3d>& bnd,
        const std::vector<std::vector<Vec3d>>& holes,
        std::vector<std::vector<Vec3d>>& outTriangles); // 每个元素是三角形（3点）

    // 是否是闭合三角网格
    static bool IsClosedSolid(const TriMesh& mesh);
};