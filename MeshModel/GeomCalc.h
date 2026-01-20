#pragma once
#include "Geometry.h"

struct MESHMODELDLL GeomCalc
{
    // 计算面片法向
    static Vec3d CompuateNormal(std::vector<Vec3d>const& ps);

    // 判断点是否在线段上
    static bool IsPointOnSegment(Vec3d const& p, Vec3d const& a, Vec3d const& b, double& paramP/*点p的参数*/);

    // 计算点到面的距离
    static double Point2PlaneDistatnce(Vec3d const& pnt, Vec3d const& o, Vec3d const& n);
    static double Point2PlaneSignDistance(Vec3d const& pnt, Vec3d const& o, Vec3d const& n);

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

    // 三角形共面裁剪
    static bool TriRegionSplit(
        std::vector<Vec3d>const& tri1, 
        std::vector<Vec3d>const& tri2, 
        std::vector<std::vector<Vec3d>>& tri1OutTri2, 
        std::vector<std::vector<Vec3d>>& tri2OutTri1, 
        std::vector<Vec3d>& tri1CollapseTri2
    );

    // 多边形的三角化（结果法向保持一致）
    static std::vector<std::vector<Vec3d>> Triangulate(std::vector<Vec3d> const& poly);

    // 包含限制边的三角化
    static std::vector<std::vector<Vec3d>> TriangulateWithConstraints(
        const std::vector<Vec3d>& outerBoundary,
        const std::vector<std::vector<Vec3d>>& intersectionPolylines
    );
};