// MeshBooleanOperation.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <vector>
#include "Geometry.h"      // 你的 Vec3d 定义
#include "GeomCalc.h"

void printPolygon(const std::vector<Vec3d>& poly, const char* name) {
    std::cout << name << " (size=" << poly.size() << "):\n";
    for (size_t i = 0; i < poly.size(); ++i) {
        std::cout << "  (" << poly[i].x << ", " << poly[i].y << ", " << poly[i].z << ")\n";
    }
}

void printMultiPolygon(const std::vector<std::vector<Vec3d>>& polys, const char* name) {
    std::cout << name << " (num polygons=" << polys.size() << "):\n";
    for (size_t i = 0; i < polys.size(); ++i) {
        std::cout << "  Polygon " << i << ":\n";
        std::vector<std::vector<Vec3d>> tris = GeomCalc::Triangulate(polys[i]);
        for (const auto& p : polys[i]) {
            std::cout << "    (" << p.x << ", " << p.y << ", " << p.z << ")\n";
        }
    }
}

int main() {
    // === 场景：两个部分重叠的共面三角形（位于 z=0 平面）===
    std::vector<Vec3d> tri1 = {
        Vec3d{0.0, 0.0, 0.0},
        Vec3d{2.0, 0.0, 0.0},
        Vec3d{1.0, 1.0, 0.0 }
    };

    std::vector<Vec3d> tri2 = {
        Vec3d{0.5, 0.0, 0.0},
        Vec3d{1.5, 0.0, 0.0 },
        Vec3d{1.0, 1.0, 0.0 }
    };

    // 输出容器
    std::vector<std::vector<Vec3d>> tri1_minus_tri2;
    std::vector<std::vector<Vec3d>> tri2_minus_tri1;
    std::vector<Vec3d> intersection;

    bool success = GeomCalc::TriRegionSplit(tri1, tri2, tri1_minus_tri2, tri2_minus_tri1, intersection);

    if (!success) {
        std::cout << "Error: Input triangles are degenerate or invalid.\n";
        return -1;
    }

    // === 打印结果 ===
    printMultiPolygon(tri1_minus_tri2, "tri1 \\ tri2");
    printMultiPolygon(tri2_minus_tri1, "tri2 \\ tri1");
    printPolygon(intersection, "tri1 ∩ tri2");

    /*
     预期结果（近似）：
     - tri1 \ tri2: 一个四边形 [(0,0,0), (1,0,0), (1,1,0), (0,2,0)]
     - tri2 \ tri1: 一个四边形 [(2,0,0), (3,0,0), (1,2,0), (1,1,0)]
     - intersection: 一个三角形 [(1,0,0), (2,0,0), (1,1,0)] 或类似（Clipper2 可能顶点顺序不同）
    */

    return 0;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
