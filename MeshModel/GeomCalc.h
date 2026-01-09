#pragma once
#include "Geometry.h"

struct MESHMODELDLL GeomCalc
{
    // 计算面片法向
    static Vec3d CompuateNormal(std::vector<Vec3d>const& ps);

    // 计算点到面的距离
    static double Point2PlaneDistatnce(Vec3d const& pnt, Vec3d const& o, Vec3d const& n);

    // 三角形共面裁剪
    static bool TriRegionSplit(
        std::vector<Vec3d>const& tri1, 
        std::vector<Vec3d>const& tri2, 
        std::vector<Vec3d>& tri1OutTri2, 
        std::vector<Vec3d>& tri2OutTri1, 
        std::vector<Vec3d>& tri1CollapseTri2
    );

    // 多边形的三角化
    static std::vector<std::vector<Vec3d>> Triangulate(std::vector<Vec3d> const& poly);
};

